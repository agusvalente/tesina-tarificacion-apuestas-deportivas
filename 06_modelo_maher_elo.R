# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 06: Modelo Maher+Elo y selección de modelo final
# =============================================================================
# Depende de: Script 04 (resultados_maher.rds)   — solo para cotejo
#             Script 05 (resultados_elo.rds)     — datos con Elo
#             resultados_lrt.rds                 — solo para cotejo de logLik
#
# -----------------------------------------------------------------------------
# DECISIONES DE DISEÑO
# -----------------------------------------------------------------------------
# (1) IDENTIFICABILIDAD. Se fijan en 0 tanto el primer ataque como la primera
#     defensa → modelo identificado, k = 50 + n_betas. Fijar solo el primer
#     ataque y dejar libres todas las defensas junto con media_global deja una
#     dirección plana exacta: sumar c a todas las defensas y restar c a
#     media_global no cambia la verosimilitud. El resultado sería un Hessiano
#     singular (sin errores estándar), k mal contado en AIC/BIC y BFGS
#     deambulando por una cresta.
#
# (2) COVARIABLES. Se testean cuatro especificaciones, no una:
#       M0 Maher puro
#       M1 + ΔElo total       (delta_elo)
#       M2 + ΔElo dinámico    (delta_elo_din)  <-- recomendada por Script 05
#       M3 + ambas
#     El ΔElo se divide por 100: en la escala original (sd ≈ 175) el
#     coeficiente resulta ~1e-4, mal condicionado para BFGS e ilegible.
#     Así β se lee como "efecto por cada 100 puntos de ΔElo".
#
# (3) CRITERIO DE SELECCIÓN. El LRT/AIC/BIC sobre train está estructuralmente
#     sesgado contra el Elo: los ~50 parámetros de equipo, estimados con los
#     mismos partidos que generaron el rating, ya absorben toda la fuerza
#     sistemática. Por eso la selección se hace por VALIDACIÓN INTERNA:
#       fit en 2021/22 + 2022/23  →  evaluar en 2023/24
#     El test (2024/25) NO participa de la selección; se reserva como
#     evaluación honesta fuera de muestra. Seleccionar con el test filtraría
#     información hacia el Script 07.
#
# (4) EVALUACIÓN PREDICTIVA. Log-loss, RPS y Brier sobre datos no usados en la
#     estimación, más diagnósticos de calibración: sesgo de localía
#     (advertido por el Script 05) y subregistro de empates (motivación de la
#     corrección de Dixon-Coles).
#
# (5) SENSIBILIDAD AL BURN-IN. Todo se estima con y sin 2021/22: en esa
#     temporada el ΔElo arranca en cero para todos los equipos y su sd es 116
#     contra 175-186 en las siguientes, por lo que la señal está atenuada.
#
# (6) VEROSIMILITUD VECTORIZADA. Un loop partido a partido haría cada
#     evaluación ~1000x más lenta. Vectorizada, el conjunto completo de
#     ajustes corre en segundos y habilita multi-arranque.
#
# -----------------------------------------------------------------------------
# ESPECIFICACIÓN
# -----------------------------------------------------------------------------
#   log(λ) = μ + ventaja_local + ataque_local     + defensa_visitante + Σ βk·Xk
#   log(μ) = μ                 + ataque_visitante + defensa_local     − Σ βk·Xk
#
# El signo es simétrico: si el local es más fuerte según el Elo, su tasa de
# goles sube y la del visitante baja en la misma magnitud.
# =============================================================================

library(tidyverse)

set.seed(904158)

# -----------------------------------------------------------------------------
# 1. CARGAR DATOS
# -----------------------------------------------------------------------------

maher <- readRDS("resultados_maher.rds")
elo   <- readRDS("resultados_elo.rds")

datos_train_elo <- elo$datos_train_elo
datos_test_elo  <- elo$datos_test_elo

equipos   <- maher$equipos          # 25 equipos presentes en train
n_equipos <- length(equipos)

cat("Equipos en el modelo (presentes en train):", n_equipos, "\n")
cat("Partidos train:", nrow(datos_train_elo),
    "| Partidos test:", nrow(datos_test_elo), "\n")

# Verificación: todos los equipos de train deben estar en 'equipos'
stopifnot(all(c(datos_train_elo$home_team, datos_train_elo$away_team) %in% equipos))

# Ipswich es el único equipo del test sin historial en train
equipos_test     <- sort(unique(c(datos_test_elo$home_team, datos_test_elo$away_team)))
equipos_nuevos   <- setdiff(equipos_test, equipos)
equipos_ausentes <- setdiff(equipos, equipos_test)

cat("Equipos nuevos en test:", paste(equipos_nuevos, collapse = ", "), "\n")
cat("Equipos de train ausentes en test:", paste(equipos_ausentes, collapse = ", "), "\n")


# -----------------------------------------------------------------------------
# 2. PREPARACIÓN DE DATOS PARA LA VEROSIMILITUD
# -----------------------------------------------------------------------------
# El ΔElo se escala dividiendo por 100 para que β sea interpretable y el
# problema esté bien condicionado.

ESCALA_ELO <- 100

# IMPORTANTE: cada ajuste usa ÚNICAMENTE los equipos presentes en su propia
# muestra de estimación. Si se usara siempre la lista completa de 25, los
# equipos ausentes de esa muestra tendrían gradiente exactamente cero: sus
# parámetros quedarían en el valor inicial, el Hessiano sería singular y k
# estaría mal contado. Ejemplos concretos en este script:
#   - muestra de validación (21/22 + 22/23): faltan Luton y Sheffield United
#   - muestra sin burn-in   (22/23 + 23/24): faltan Norwich y Watford
# Los equipos faltantes al momento de predecir se imputan igual que Ipswich.

preparar <- function(d, covariables, eq) {
  list(
    il   = match(d$home_team, eq),
    iv   = match(d$away_team, eq),
    gl   = d$fthg,
    gv   = d$ftag,
    ftr  = d$ftr,
    X    = if (length(covariables) == 0) NULL else
             as.matrix(d[, covariables, drop = FALSE]) / ESCALA_ELO,
    n    = nrow(d)
  )
}

# Definición de las especificaciones a comparar
ESPECIFICACIONES <- list(
  M0 = list(nombre = "Maher puro",        covs = character(0)),
  M1 = list(nombre = "Maher + ΔElo",      covs = "delta_elo"),
  M2 = list(nombre = "Maher + ΔElo din.", covs = "delta_elo_din"),
  M3 = list(nombre = "Maher + ambas",     covs = c("delta_elo_ini", "delta_elo_din"))
)

# Nota sobre M3: se usan delta_elo_ini y delta_elo_din en lugar de
# delta_elo y delta_elo_din porque delta_elo = ini + din exactamente, y meter
# las tres sería colinealidad perfecta. La descomposición (ini, din) es la
# parametrización equivalente y no degenerada.


# -----------------------------------------------------------------------------
# 3. VEROSIMILITUD VECTORIZADA
# -----------------------------------------------------------------------------
# Disposición de parámetros:
#   [1]                media_global
#   [2]                ventaja_local
#   [3 : 2+nb]         betas
#   [3+nb : 1+nb+ne]   ataques  de los equipos 2..ne   (el 1 queda fijo en 0)
#   [2+nb+ne : nb+2ne] defensas de los equipos 2..ne   (la 1 queda fija en 0)
#
# Ambas restricciones son necesarias: fijar solo los ataques deja la dirección
# plana (media_global, defensas) descrita en el encabezado.

desarmar_params <- function(params, nb, ne) {
  off <- 2 + nb
  list(
    media_global  = params[1],
    ventaja_local = params[2],
    beta          = if (nb > 0) params[3:(2 + nb)] else numeric(0),
    ataques       = c(0, params[(off + 1):(off + ne - 1)]),
    defensas      = c(0, params[(off + ne):(off + 2 * (ne - 1))])
  )
}

neg_loglik <- function(params, dat, nb, ne) {
  p <- desarmar_params(params, nb, ne)

  eta <- if (nb > 0) as.numeric(dat$X %*% p$beta) else 0

  log_lambda <- p$media_global + p$ventaja_local +
                p$ataques[dat$il] + p$defensas[dat$iv] + eta
  log_mu     <- p$media_global +
                p$ataques[dat$iv] + p$defensas[dat$il] - eta

  # Cota de seguridad: evita Inf si el optimizador explora valores extremos
  if (any(!is.finite(log_lambda)) || any(!is.finite(log_mu))) return(1e10)
  if (max(log_lambda, log_mu) > 5) return(1e10)

  -sum(dpois(dat$gl, exp(log_lambda), log = TRUE) +
       dpois(dat$gv, exp(log_mu),     log = TRUE))
}


# -----------------------------------------------------------------------------
# 4. FUNCIÓN DE AJUSTE (con multi-arranque)
# -----------------------------------------------------------------------------

ajustar <- function(datos, covs, n_arranques = 3, verbose = FALSE) {

  # Equipos efectivamente presentes en ESTA muestra (ver nota en la sección 2)
  eq  <- sort(unique(c(datos$home_team, datos$away_team)))
  ne  <- length(eq)
  nb  <- length(covs)
  dat <- preparar(datos, covs, eq)
  np  <- nb + 2 * ne

  mejor <- NULL

  for (s in seq_len(n_arranques)) {
    init <- c(
      log(mean(datos$fthg)),            # media global
      0.15,                             # ventaja local
      rep(0, nb),                       # betas
      rep(0, 2 * (ne - 1))              # ataques y defensas
    )
    # Arranques 2 y 3: perturbación aleatoria para chequear óptimos locales
    if (s > 1) init <- init + rnorm(np, 0, 0.05)

    fit <- try(optim(
      par     = init,
      fn      = neg_loglik,
      dat     = dat, nb = nb, ne = ne,
      method  = "BFGS",
      control = list(maxit = 10000, reltol = 1e-12),
      hessian = TRUE
    ), silent = TRUE)

    if (inherits(fit, "try-error")) next
    if (is.null(mejor) || fit$value < mejor$value) mejor <- fit
  }

  stopifnot(!is.null(mejor))

  p <- desarmar_params(mejor$par, nb, ne)

  # Errores estándar de los betas (el Hessiano es invertible, ver (1))
  se_beta <- rep(NA_real_, nb)
  if (nb > 0) {
    vcov <- try(solve(mejor$hessian), silent = TRUE)
    if (!inherits(vcov, "try-error")) {
      d <- diag(vcov)[3:(2 + nb)]
      se_beta <- ifelse(d > 0, sqrt(d), NA_real_)
    }
  }

  # Reparametrización a suma cero (solo para reportar; no altera λ ni μ)
  a_c <- mean(p$ataques); d_c <- mean(p$defensas)
  ataques_c  <- p$ataques  - a_c
  defensas_c <- p$defensas - d_c
  media_c    <- p$media_global + a_c + d_c

  names(ataques_c) <- names(defensas_c) <- eq

  list(
    covs          = covs,
    nb            = nb,
    k             = np,
    equipos       = eq,
    n_equipos     = ne,
    loglik        = -mejor$value,
    convergencia  = mejor$convergence,
    par           = mejor$par,
    beta          = p$beta,
    se_beta       = se_beta,
    media_global  = media_c,
    ventaja_local = p$ventaja_local,
    ataques       = ataques_c,
    defensas      = defensas_c,
    n_obs         = nrow(datos)
  )
}


# -----------------------------------------------------------------------------
# 5. PROBABILIDADES 1X2 (vectorizado)
# -----------------------------------------------------------------------------

MAX_GOLES <- 25   # mismo valor que el Script 04

probs_1x2 <- function(lambda, mu, max_goles = MAX_GOLES) {
  g <- 0:max_goles

  # OJO CON EL ORDEN DE ARGUMENTOS: outer(X, Y, f) evalúa f(X[i], Y[j]).
  # Escribir outer(lambda, g, dpois) llama dpois(x = lambda, lambda = g), es
  # decir pide la probabilidad de observar un número NO ENTERO de goles bajo
  # una Poisson de media 0,1,2,... dpois devuelve 0 con warning, toda la
  # matriz queda nula y las probabilidades salen 0/0 = NaN.
  # La lambda anónima fuerza el orden correcto: dpois(goles, tasa).
  PL <- outer(lambda, g, \(tasa, goles) dpois(goles, tasa))  # n x (G+1)
  PV <- outer(mu,     g, \(tasa, goles) dpois(goles, tasa))

  cumPL <- t(apply(PL, 1, cumsum))
  cumPV <- t(apply(PV, 1, cumsum))

  # P(H) = Σ_i P(local = i) · P(visitante < i)
  menorPV <- cbind(0, cumPV[, -ncol(cumPV), drop = FALSE])
  menorPL <- cbind(0, cumPL[, -ncol(cumPL), drop = FALSE])

  pH <- rowSums(PL * menorPV)
  pD <- rowSums(PL * PV)
  pA <- rowSums(PV * menorPL)

  tot <- pH + pD + pA                      # renormaliza la cola truncada
  P <- cbind(prob_H = pH / tot, prob_D = pD / tot, prob_A = pA / tot)

  # Red de seguridad: cualquier NaN o suma distinta de 1 indica un problema
  # aguas arriba (tasas explotadas, parámetros NA). Mejor fallar acá.
  stopifnot(all(is.finite(P)), max(abs(rowSums(P) - 1)) < 1e-8)
  P
}

# Predicción para un conjunto de partidos usando un modelo ajustado.
# Un equipo puede faltar en el ajuste por dos motivos equivalentes:
#   - nunca jugó en la muestra de estimación (Ipswich frente a train)
#   - ascendió después del corte (Luton y Sheffield United frente a la
#     muestra de validación; Norwich y Watford frente a la muestra sin
#     burn-in)
# En todos los casos se aplica el criterio del Script 04: media ponderada de
# los parámetros ajustados, con peso 3 sobre los equipos que estaban en la
# muestra de estimación pero no en la de predicción (los "descendidos").
PESO_DESCENDIDOS <- 3

predecir <- function(fit, datos, verbose = FALSE) {

  eq_fit   <- fit$equipos
  eq_pred  <- sort(unique(c(datos$home_team, datos$away_team)))
  nuevos   <- setdiff(eq_pred, eq_fit)
  ausentes <- setdiff(eq_fit,  eq_pred)

  ataques  <- fit$ataques
  defensas <- fit$defensas

  if (length(nuevos) > 0) {
    pesos <- ifelse(eq_fit %in% ausentes, PESO_DESCENDIDOS, 1)
    a_imp <- weighted.mean(fit$ataques,  pesos)
    d_imp <- weighted.mean(fit$defensas, pesos)
    for (e in nuevos) {
      ataques[e]  <- a_imp
      defensas[e] <- d_imp
    }
    if (verbose) {
      cat("   Equipos imputados:", paste(nuevos, collapse = ", "),
          "| ataque =", round(a_imp, 4), "| defensa =", round(d_imp, 4), "\n")
    }
  }

  eta <- 0
  if (fit$nb > 0) {
    X   <- as.matrix(datos[, fit$covs, drop = FALSE]) / ESCALA_ELO
    eta <- as.numeric(X %*% fit$beta)
  }

  a_l <- ataques[datos$home_team];  d_v <- defensas[datos$away_team]
  a_v <- ataques[datos$away_team];  d_l <- defensas[datos$home_team]
  stopifnot(!any(is.na(c(a_l, d_v, a_v, d_l))))

  lambda <- as.numeric(exp(fit$media_global + fit$ventaja_local + a_l + d_v + eta))
  mu     <- as.numeric(exp(fit$media_global + a_v + d_l - eta))

  P <- probs_1x2(lambda, mu)

  datos |>
    mutate(lambda = lambda,
           mu     = mu,
           prob_H = P[, "prob_H"],
           prob_D = P[, "prob_D"],
           prob_A = P[, "prob_A"])
}


# -----------------------------------------------------------------------------
# 6. MÉTRICAS PREDICTIVAS
# -----------------------------------------------------------------------------
# log-loss : verosimilitud predictiva multinomial (menor = mejor)
# RPS      : Ranked Probability Score, penaliza según la distancia ordinal
#            entre el resultado previsto y el observado (H > D > A).
#            Es la métrica estándar en pronóstico de fútbol.
# Brier    : error cuadrático multiclase

# Pérdida PARTIDO A PARTIDO: necesaria para los contrastes pareados.
perdidas <- function(pred) {
  pH <- pred$prob_H; pD <- pred$prob_D; pA <- pred$prob_A
  oH <- as.numeric(pred$ftr == "H")
  oD <- as.numeric(pred$ftr == "D")
  oA <- as.numeric(pred$ftr == "A")

  p_obs <- pH * oH + pD * oD + pA * oA

  tibble(
    ll    = -log(pmax(p_obs, 1e-12)),
    rps_i = ((pH - oH)^2 + ((pH + pD) - (oH + oD))^2) / 2,
    br_i  = (pH - oH)^2 + (pD - oD)^2 + (pA - oA)^2,
    ok    = as.numeric(max.col(cbind(pH, pD, pA)) == max.col(cbind(oH, oD, oA)))
  )
}

metricas <- function(pred) {
  L <- perdidas(pred)
  tibble(
    n        = nrow(pred),
    log_loss = mean(L$ll),
    rps      = mean(L$rps_i),
    brier    = mean(L$br_i),
    acierto  = mean(L$ok)
  )
}

# -----------------------------------------------------------------------------
# CONTRASTE PAREADO (Diebold-Mariano)
# -----------------------------------------------------------------------------
# Una diferencia de log-loss de 0.007 sobre 380 partidos puede ser ruido puro.
# El contraste correcto opera sobre las diferencias PARTIDO A PARTIDO de la
# pérdida: d_i = pérdida(modelo) − pérdida(M0). Si E[d] < 0 de forma
# significativa, el modelo predice mejor. Al ser pareado, elimina la varianza
# común de los partidos y tiene mucha más potencia que comparar dos medias
# independientes.
#
# Nota: los modelos están anidados y estimados sobre los mismos datos, así que
# el DM estándar es conservador (Clark & West). Se usa como criterio de corte
# práctico, no como test exacto.

comparar_pareado <- function(fits, datos_eval, etiqueta, referencia = "M0") {
  L <- map(fits, \(f) perdidas(predecir(f, datos_eval)))

  imap_dfr(L, \(li, id) {
    if (id == referencia) {
      return(tibble(muestra = etiqueta, modelo = id,
                    dif_log_loss = 0, ic_inf = NA_real_, ic_sup = NA_real_,
                    p_DM = NA_real_, dif_rps = 0, p_DM_rps = NA_real_))
    }
    d_ll  <- li$ll    - L[[referencia]]$ll
    d_rps <- li$rps_i - L[[referencia]]$rps_i
    t_ll  <- t.test(d_ll)
    t_rps <- t.test(d_rps)
    tibble(
      muestra      = etiqueta,
      modelo       = id,
      dif_log_loss = mean(d_ll),
      ic_inf       = t_ll$conf.int[1],
      ic_sup       = t_ll$conf.int[2],
      p_DM         = t_ll$p.value,
      dif_rps      = mean(d_rps),
      p_DM_rps     = t_rps$p.value
    )
  })
}


# -----------------------------------------------------------------------------
# 7. MUESTRAS DE ESTIMACIÓN
# -----------------------------------------------------------------------------

temporadas_train <- sort(unique(datos_train_elo$temporada))
primera_temp     <- temporadas_train[1]
ultima_train     <- temporadas_train[length(temporadas_train)]

# (a) Train completo
train_full <- datos_train_elo

# (b) Train sin la primera temporada (burn-in del Elo)
train_sin_burnin <- datos_train_elo |> filter(temporada != primera_temp)

# (c) Validación interna leakage-free: se ajusta con todo menos la última
#     temporada de train y se evalúa sobre ella.
val_fit <- datos_train_elo |> filter(temporada != ultima_train)
val_eval <- datos_train_elo |> filter(temporada == ultima_train)

cat("\n--- Muestras ---\n")
cat("Train completo      :", nrow(train_full), "partidos\n")
cat("Train sin", primera_temp, ":", nrow(train_sin_burnin), "partidos\n")
cat("Validación: fit =", nrow(val_fit), "| eval =", nrow(val_eval),
    "(", ultima_train, ")\n")

# Diagnóstico de cobertura de equipos por muestra. Los que aparecen como
# "a imputar" son equipos que juegan en la muestra de evaluación pero no
# estaban en la de estimación: reciben el tratamiento de Ipswich.
cobertura <- function(fit_d, eval_d, etiqueta) {
  eqf <- sort(unique(c(fit_d$home_team,  fit_d$away_team)))
  eqe <- sort(unique(c(eval_d$home_team, eval_d$away_team)))
  cat(sprintf("%-28s | equipos estimados: %2d | a imputar: %s\n",
              etiqueta, length(eqf),
              if (length(setdiff(eqe, eqf)) == 0) "ninguno"
              else paste(setdiff(eqe, eqf), collapse = ", ")))
}

cat("\n--- Cobertura de equipos ---\n")
cobertura(train_full,       datos_test_elo, "train completo → test")
cobertura(train_sin_burnin, datos_test_elo, paste("sin", primera_temp, "→ test"))
cobertura(val_fit,          val_eval,       "validación → 2023/24")


# -----------------------------------------------------------------------------
# 8. AJUSTE IN-SAMPLE Y COMPARACIÓN LRT / AIC / BIC
# -----------------------------------------------------------------------------
# Se reporta pero NO decide: es la métrica sesgada contra el Elo.
# Nota: los modelos se reestiman todos aquí con la MISMA parametrización y el
# MISMO optimizador, de modo que la diferencia de log-verosimilitud no puede
# deberse a artefactos entre corridas distintas.

ajustar_todas <- function(datos, etiqueta) {
  cat("\nAjustando especificaciones sobre", etiqueta,
      "(", nrow(datos), "partidos )...\n")
  fits <- map(ESPECIFICACIONES, \(e) ajustar(datos, e$covs))
  names(fits) <- names(ESPECIFICACIONES)
  fits
}

tabla_insample <- function(fits, datos, etiqueta) {
  n  <- nrow(datos)
  m0 <- fits$M0

  imap_dfr(fits, \(f, id) {
    LR <- 2 * (f$loglik - m0$loglik)
    df <- f$k - m0$k
    tibble(
      muestra   = etiqueta,
      modelo    = id,
      nombre    = ESPECIFICACIONES[[id]]$nombre,
      k         = f$k,
      loglik    = f$loglik,
      AIC       = -2 * f$loglik + 2 * f$k,
      BIC       = -2 * f$loglik + f$k * log(n),
      LR        = if (df > 0) LR else NA_real_,
      gl        = if (df > 0) df else NA_integer_,
      p_LRT     = if (df > 0) pchisq(LR, df, lower.tail = FALSE) else NA_real_,
      converge  = f$convergencia == 0
    )
  })
}

fits_full   <- ajustar_todas(train_full,       "train completo")
fits_sinbi  <- ajustar_todas(train_sin_burnin, paste("train sin", primera_temp))

insample <- bind_rows(
  tabla_insample(fits_full,  train_full,       "train completo"),
  tabla_insample(fits_sinbi, train_sin_burnin, paste("sin", primera_temp))
)

cat("\n=== AJUSTE IN-SAMPLE (referencia, no decide) ===\n")
insample |>
  mutate(across(c(loglik, AIC, BIC, LR), \(x) round(x, 2)),
         p_LRT = round(p_LRT, 4)) |>
  print(width = Inf)

cat("\n>> El LRT sobre train está sesgado CONTRA el Elo: los", 2 * (n_equipos - 1),
    "parámetros\n")
cat("   de equipo se estiman con los mismos partidos que generaron el rating y\n")
cat("   ya absorben la fuerza sistemática. La decisión se toma en la sección 9.\n")

# --- Coeficientes beta con error estándar -----------------------------------
cat("\n=== COEFICIENTES β (train completo) ===\n")
imap_dfr(fits_full, \(f, id) {
  if (f$nb == 0) return(NULL)
  tibble(
    modelo    = id,
    covariable = f$covs,
    beta      = f$beta,
    se        = f$se_beta,
    z         = f$beta / f$se_beta,
    p         = 2 * pnorm(-abs(f$beta / f$se_beta))
  )
}) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print(width = Inf)

cat("\n>> β está expresado por cada", ESCALA_ELO, "puntos de ΔElo.\n")

# --- Cotejo con los scripts previos -----------------------------------------
cat("\n--- Cotejo con Script 04 / resultados_lrt.rds ---\n")
if (file.exists("resultados_lrt.rds")) {
  loglik_prev <- readRDS("resultados_lrt.rds")$loglik_maher
  cat("logLik Maher (script previo):", round(loglik_prev, 3), "\n")
  cat("logLik Maher (reestimado)   :", round(fits_full$M0$loglik, 3), "\n")
  cat("Diferencia                  :",
      round(fits_full$M0$loglik - loglik_prev, 4), "\n")
  if (abs(fits_full$M0$loglik - loglik_prev) > 0.05) {
    cat(">> AVISO: discrepancia no trivial. Alguna de las dos optimizaciones no\n")
    cat("   alcanzó el máximo. Se usa la de este script (multi-arranque,\n")
    cat("   parametrización identificada).\n")
  } else {
    cat(">> Coinciden: la reestimación reproduce el ajuste previo.\n")
  }
} else {
  cat("resultados_lrt.rds no encontrado; se omite el cotejo.\n")
}


# -----------------------------------------------------------------------------
# 9. SELECCIÓN POR VALIDACIÓN INTERNA (criterio decisorio)
# -----------------------------------------------------------------------------
# Se ajusta con 2021/22 + 2022/23 y se evalúa sobre 2023/24. Es la única
# comparación que reproduce el uso real del modelo (predecir una temporada no
# vista) sin tocar el test. La 2024/25 queda intacta para la evaluación final.

cat("\n=== VALIDACIÓN INTERNA: fit hasta", ultima_train, "excluida, eval sobre",
    ultima_train, "===\n")

fits_val <- map(ESPECIFICACIONES, \(e) ajustar(val_fit, e$covs))
names(fits_val) <- names(ESPECIFICACIONES)

cat("Equipos estimados en la muestra de validación:",
    fits_val$M0$n_equipos, "\n")

validacion <- imap_dfr(fits_val, \(f, id) {
  pred <- predecir(f, val_eval, verbose = (id == "M0"))
  bind_cols(
    tibble(modelo = id, nombre = ESPECIFICACIONES[[id]]$nombre),
    metricas(pred)
  )
})

stopifnot(all(is.finite(validacion$log_loss)))

base_ll <- validacion$log_loss[validacion$modelo == "M0"]

validacion <- validacion |>
  mutate(delta_log_loss = log_loss - base_ll,
         mejora_pct     = 100 * (base_ll - log_loss) / base_ll)

validacion |>
  mutate(across(where(is.numeric), \(x) round(x, 5))) |>
  print(width = Inf)

# --- Contraste pareado: ¿la diferencia es distinguible de cero? -------------
# Sin esto, "M1 dio 0.938 y M0 dio 0.945" es un número sin interpretación:
# no se sabe si la brecha excede el ruido muestral de 380 partidos.

cat("\n--- Contraste pareado contra M0 (validación) ---\n")
dm_val <- comparar_pareado(fits_val, val_eval, "validación")
dm_val |>
  mutate(across(where(is.numeric), \(x) round(x, 5))) |>
  print(width = Inf)

cat("\n>> dif_log_loss negativo = el modelo predice mejor que M0.\n")
cat(">> p_DM es el p-valor del contraste pareado sobre las diferencias por\n")
cat("   partido. Si el intervalo de confianza cruza el cero, la mejora no es\n")
cat("   distinguible del ruido y no debe reportarse como un hallazgo.\n")

# --- Regla de selección -----------------------------------------------------
# Se recorren las especificaciones ordenadas por log-loss de validación y se
# adopta la primera que cumpla LOS TRES requisitos:
#
#   (a) menor log-loss que M0
#   (b) contraste pareado que respalde la mejora (α = 0.10, laxo porque el DM
#       sobre modelos anidados es conservador)
#   (c) ESTABILIDAD DE SIGNO de β entre las dos muestras de estimación
#
# El requisito (c) es una exigencia estándar de estabilidad paramétrica, no una
# regla pre-registrada: se incorporó en el curso del análisis, al observar la
# inversión de signo. Opera enteramente sobre las muestras de estimación, sin
# consultar en ningún momento el test, de modo que no constituye selección
# oportunista sobre el resultado. Un coeficiente que se invierte de signo al
# cambiar la muestra no está identificado por los datos: está absorbiendo
# ruido por colinealidad con los parámetros de equipo. Un modelo así puede
# predecir mejor por casualidad en una evaluación puntual, pero su parámetro
# no admite interpretación y no es defendible. La estabilidad del parámetro
# prima sobre una diferencia marginal de log-loss.
#
# Contexto de por qué esto importa acá: en 2021/22 el delta_elo_ini es cero
# para todos los equipos, de modo que delta_elo es puro componente dinámico y
# NO es colineal con los parámetros de Maher. Esas 380 observaciones son las
# que aportan la variación identificante. Al excluirlas queda solo data donde
# el ΔElo es mayoritariamente nivel estático, casi indistinguible de los
# parámetros de equipo, y el coeficiente se vuelve inestimable.

ALFA_DM <- 0.10

signo_estable <- function(id) {
  if (length(ESPECIFICACIONES[[id]]$covs) == 0) return(TRUE)
  b_con <- fits_full[[id]]$beta
  b_sin <- fits_sinbi[[id]]$beta
  all(sign(b_con) == sign(b_sin))
}

orden_candidatos <- validacion |> arrange(log_loss) |> pull(modelo)

cat("\n--- Evaluación de candidatos (orden por log-loss de validación) ---\n")

modelo_sel <- "M0"
for (id in orden_candidatos) {

  if (id == "M0") {
    cat("M0: referencia. Ningún candidato previo cumplió los tres requisitos.\n")
    modelo_sel <- "M0"
    break
  }

  p_id <- dm_val$p_DM[dm_val$modelo == id]
  d_id <- dm_val$dif_log_loss[dm_val$modelo == id]
  est  <- signo_estable(id)

  b_con <- fits_full[[id]]$beta
  b_sin <- fits_sinbi[[id]]$beta

  cat(sprintf("\n%s (%s)\n", id, ESPECIFICACIONES[[id]]$nombre))
  cat("  (a) dif log-loss vs M0:", round(d_id, 5),
      if (d_id < 0) " ✓" else " ✗", "\n")
  cat("  (b) p contraste pareado:", round(p_id, 5),
      if (!is.na(p_id) && p_id < ALFA_DM) " ✓" else " ✗", "\n")
  cat("  (c) signo de β estable :",
      paste(round(b_con, 4), "vs", round(b_sin, 4), collapse = " | "),
      if (est) " ✓" else " ✗ SE INVIERTE", "\n")

  if (d_id < 0 && !is.na(p_id) && p_id < ALFA_DM && est) {
    modelo_sel <- id
    cat("  → CUMPLE los tres requisitos. Seleccionado.\n")
    break
  } else {
    cat("  → descartado.\n")
  }
}

cat("\n=== MODELO SELECCIONADO:", modelo_sel, "—",
    ESPECIFICACIONES[[modelo_sel]]$nombre, "===\n")
cat(">> Criterio: menor log-loss en validación interna, respaldado por\n")
cat("   contraste pareado (α =", ALFA_DM, ") y con signo de β estable entre\n")
cat("   muestras de estimación.\n")
cat(">> El test 2024/25 NO participó de esta decisión.\n")


# -----------------------------------------------------------------------------
# 10. SENSIBILIDAD AL BURN-IN EN LA MUESTRA DE ESTIMACIÓN
# -----------------------------------------------------------------------------
# ¿Cambia el ranking de modelos si se estima sin 2021/22? El ΔElo de esa
# temporada arranca en cero para todos (sd 116 vs 175-186), señal atenuada.

cat("\n=== SENSIBILIDAD: estimación con vs sin", primera_temp, "===\n")

sensibilidad <- bind_rows(
  imap_dfr(fits_full, \(f, id) bind_cols(
    tibble(muestra = "train completo", modelo = id), metricas(predecir(f, datos_test_elo))
  )),
  imap_dfr(fits_sinbi, \(f, id) bind_cols(
    tibble(muestra = paste("sin", primera_temp), modelo = id),
    metricas(predecir(f, datos_test_elo))
  ))
)

sensibilidad |>
  select(muestra, modelo, log_loss, rps, brier) |>
  pivot_wider(names_from = muestra, values_from = c(log_loss, rps, brier)) |>
  mutate(across(where(is.numeric), \(x) round(x, 5))) |>
  print(width = Inf)

# Alerta automática de fragilidad: ¿el modelo elegido conserva su ventaja
# cuando se estima sin la temporada de burn-in?
ll_wide <- sensibilidad |>
  select(muestra, modelo, log_loss) |>
  pivot_wider(names_from = muestra, values_from = log_loss)

col_con <- "train completo"
col_sin <- paste("sin", primera_temp)

rank_con <- ll_wide$modelo[order(ll_wide[[col_con]])]
rank_sin <- ll_wide$modelo[order(ll_wide[[col_sin]])]

cat("\nRanking con", primera_temp, ":", paste(rank_con, collapse = " > "), "\n")
cat("Ranking sin", primera_temp, ":", paste(rank_sin, collapse = " > "), "\n")

if (modelo_sel != "M0") {
  d_con <- ll_wide[[col_con]][ll_wide$modelo == modelo_sel] -
           ll_wide[[col_con]][ll_wide$modelo == "M0"]
  d_sin <- ll_wide[[col_sin]][ll_wide$modelo == modelo_sel] -
           ll_wide[[col_sin]][ll_wide$modelo == "M0"]
  cat("\nVentaja de", modelo_sel, "sobre M0 en test:\n")
  cat("  estimado con", primera_temp, ":", round(d_con, 5), "\n")
  cat("  estimado sin", primera_temp, ":", round(d_sin, 5), "\n")

  # β estimado en cada muestra. Comparar SIGNO antes que magnitud: una
  # inversión de signo indica que el coeficiente no está identificado.
  b_con <- fits_full[[modelo_sel]]$beta
  b_sin <- fits_sinbi[[modelo_sel]]$beta
  s_con <- fits_full[[modelo_sel]]$se_beta
  s_sin <- fits_sinbi[[modelo_sel]]$se_beta

  cat("\nβ estimado según la muestra:\n")
  cat("  con", primera_temp, "(", fits_full[[modelo_sel]]$n_obs, "partidos ):",
      paste(round(b_con, 4), "±", round(s_con, 4), collapse = " | "), "\n")
  cat("  sin", primera_temp, "(", fits_sinbi[[modelo_sel]]$n_obs, "partidos ):",
      paste(round(b_sin, 4), "±", round(s_sin, 4), collapse = " | "), "\n")

  if (!all(sign(b_con) == sign(b_sin))) {
    cat("\n>> El signo de β SE INVIERTE entre muestras. Ver sección 9: este\n")
    cat("   criterio ya se aplicó en la selección.\n")
  } else if (d_con < 0 && d_sin > 0) {
    cat("\n>> El signo de β se mantiene pero su magnitud varía lo suficiente\n")
    cat("   como para revertir la ventaja fuera de muestra. Reportar la\n")
    cat("   imprecisión del tamaño del efecto como limitación.\n")
  } else if (d_sin < 0) {
    cat("\n>> La ventaja se mantiene en ambas muestras: conclusión robusta\n")
    cat("   tanto en signo como en magnitud.\n")
  }
} else {
  cat("\n>> Modelo seleccionado M0 (Maher puro).\n")
  cat(">> Diagnóstico de por qué ninguna especificación con Elo prosperó:\n")
  for (id in c("M1", "M2", "M3")) {
    b1 <- fits_full[[id]]$beta;  b2 <- fits_sinbi[[id]]$beta
    s1 <- fits_full[[id]]$se_beta
    cat(sprintf("   %s: β(train completo) = %s ± %s | β(sin %s) = %s%s\n",
                id,
                paste(round(b1, 4), collapse = ", "),
                paste(round(s1, 4), collapse = ", "),
                primera_temp,
                paste(round(b2, 4), collapse = ", "),
                if (!all(sign(b1) == sign(b2))) "   <-- signo invertido" else ""))
  }
  cat("   El ΔElo es en gran medida colineal con los parámetros de ataque y\n")
  cat("   defensa de Maher (ver Script 05, sección 8.4), de modo que su\n")
  cat("   coeficiente no queda identificado de forma estable por los datos.\n")
}


# -----------------------------------------------------------------------------
# 11. MODELO FINAL Y EVALUACIÓN FUERA DE MUESTRA
# -----------------------------------------------------------------------------
# El modelo seleccionado se reestima sobre TODO el train y se evalúa sobre la
# 2024/25. Estos números se REPORTAN, no se usan para elegir nada.

fit_final <- fits_full[[modelo_sel]]

resultados_fin <- predecir(fit_final, datos_test_elo)

cat("\n=== EVALUACIÓN FUERA DE MUESTRA (2024/25) ===\n")
cat("Modelo final:", ESPECIFICACIONES[[modelo_sel]]$nombre, "\n\n")

oos <- imap_dfr(fits_full, \(f, id) bind_cols(
  tibble(modelo = id, nombre = ESPECIFICACIONES[[id]]$nombre),
  metricas(predecir(f, datos_test_elo))
))

oos |>
  mutate(seleccionado = if_else(modelo == modelo_sel, "<--", ""),
         across(where(is.numeric), \(x) round(x, 5))) |>
  print(width = Inf)

cat("\n--- Contraste pareado contra M0 (test 2024/25) ---\n")
dm_oos <- comparar_pareado(fits_full, datos_test_elo, "test 2024/25")
dm_oos |>
  mutate(across(where(is.numeric), \(x) round(x, 5))) |>
  print(width = Inf)

cat("\n>> Este contraste se REPORTA, no decide nada: la selección ya se hizo\n")
cat("   en la sección 9 con la validación interna.\n")
cat(">> Comparar el ranking de esta tabla con el de la validación. Si difieren,\n")
cat("   es porque 380 partidos dan poca resolución y la diferencia entre\n")
cat("   modelos es del orden del ruido muestral. Decirlo explícitamente.\n")


# -----------------------------------------------------------------------------
# 12. DIAGNÓSTICOS DE CALIBRACIÓN DEL MODELO FINAL
# -----------------------------------------------------------------------------

cat("\n=== 12.1 Calibración agregada por resultado ===\n")

calib_agregada <- resultados_fin |>
  summarise(
    pred_H = mean(prob_H), obs_H = mean(ftr == "H"),
    pred_D = mean(prob_D), obs_D = mean(ftr == "D"),
    pred_A = mean(prob_A), obs_A = mean(ftr == "A")
  ) |>
  pivot_longer(everything(),
               names_to = c(".value", "resultado"),
               names_sep = "_") |>
  mutate(brecha = pred - obs,
         across(where(is.numeric), \(x) round(x, 4)))

print(calib_agregada)

# Contraste de calibración: cada partido tiene su propia probabilidad p_i, de
# modo que el número de ocurrencias sigue una Poisson-binomial, no una
# binomial con p constante. Usar binom.test(p = media) sobrestima la varianza
# (Var Poisson-binomial = Σ p_i(1−p_i) ≤ n·p̄(1−p̄)) y vuelve el test
# innecesariamente conservador. Se usa el z de Poisson-binomial.

cat("\nContraste de calibración (z de Poisson-binomial):\n")
calib_test <- map_dfr(c("H", "D", "A"), \(r) {
  p_i   <- resultados_fin[[paste0("prob_", r)]]
  obs_n <- sum(resultados_fin$ftr == r)
  esp_n <- sum(p_i)
  sd_n  <- sqrt(sum(p_i * (1 - p_i)))
  z     <- (obs_n - esp_n) / sd_n
  tibble(resultado = r,
         observado = obs_n,
         esperado  = round(esp_n, 1),
         sd        = round(sd_n, 2),
         z         = round(z, 3),
         p         = round(2 * pnorm(-abs(z)), 4))
})
print(calib_test)

cat("\n>> SESGO DE LOCALÍA (advertido por el Script 05): si pred_H > obs_H de\n")
cat("   forma significativa, el modelo sobreestima al local en todo el test.\n")
cat("   En el Script 07 eso se traduce en señales de 'valor' espurias en\n")
cat("   apuestas al local. Reportar el rendimiento por lado apostado.\n")

cat("\n>> SUBREGISTRO DE EMPATES: si pred_D < obs_D, es la limitación conocida\n")
cat("   del Poisson independiente. Es lo que corrige el término τ de\n")
cat("   Dixon-Coles (1997), que infla los resultados de bajo marcador.\n")


cat("\n=== 12.2 Calibración por deciles de P(H) ===\n")

resultados_fin |>
  mutate(bin = cut(prob_H, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(n = n(),
            esperado  = round(mean(prob_H), 3),
            observado = round(mean(ftr == "H"), 3),
            .groups = "drop") |>
  print()


cat("\n=== 12.3 Marcadores de bajo score (motivación Dixon-Coles) ===\n")
# El Poisson independiente subestima 0-0 y 1-1 y sobreestima 1-0 y 0-1.
# Se compara la frecuencia observada con la esperada bajo el modelo.

marcadores <- tibble(gl = c(0, 1, 0, 1), gv = c(0, 0, 1, 1)) |>
  mutate(
    marcador  = paste0(gl, "-", gv),
    observado = map2_dbl(gl, gv, \(a, b)
                  mean(resultados_fin$fthg == a & resultados_fin$ftag == b)),
    esperado  = map2_dbl(gl, gv, \(a, b)
                  mean(dpois(a, resultados_fin$lambda) *
                       dpois(b, resultados_fin$mu)))
  ) |>
  mutate(brecha = round(observado - esperado, 4),
         across(c(observado, esperado), \(x) round(x, 4))) |>
  select(marcador, observado, esperado, brecha)

print(marcadores)

cat("\n>> Patrón típico que justifica Dixon-Coles: 0-0 y 1-1 subestimados,\n")
cat("   1-0 y 0-1 sobreestimados. Si aparece, documentarlo como limitación\n")
cat("   (o implementar la corrección τ como extensión).\n")


# -----------------------------------------------------------------------------
# 13. GRÁFICOS
# -----------------------------------------------------------------------------

# --- 13.1 Calibración del modelo final --------------------------------------
g_calib <- resultados_fin |>
  mutate(bin = cut(prob_H, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(n = n(), esperado = mean(prob_H),
            observado = mean(ftr == "H"), .groups = "drop") |>
  ggplot(aes(x = esperado, y = observado, size = n)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(alpha = 0.8, color = "#2980b9") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title    = paste0("Calibración de P(victoria local) — ",
                         ESPECIFICACIONES[[modelo_sel]]$nombre),
       subtitle = "Temporada 2024/25 | Puntos bajo la diagonal = sobreestimación del local",
       x = "Probabilidad predicha", y = "Frecuencia observada", size = "n") +
  theme_minimal()

print(g_calib)


# --- 13.2 Cuánto mueve el Elo las probabilidades ----------------------------
# Compara M0 contra el modelo seleccionado. Muestra la magnitud práctica del
# ajuste que introduce el Elo, que es distinta de su significatividad
# estadística: un efecto puede ser detectable y aun así mover poco las
# probabilidades.

pred_M0  <- predecir(fits_full$M0, datos_test_elo)
pred_sel <- resultados_fin
comparable <- modelo_sel != "M0"

if (comparable) {
  dif_p <- pred_sel$prob_H - pred_M0$prob_H

  g_comparacion <- tibble(p_M0 = pred_M0$prob_H, p_sel = pred_sel$prob_H,
                          dif = dif_p) |>
    ggplot(aes(x = p_M0, y = p_sel, color = abs(dif))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#e74c3c") +
    geom_point(alpha = 0.6, size = 1.6) +
    scale_color_viridis_c(option = "viridis", direction = -1) +
    coord_equal() +
    labs(title    = paste0("P(victoria local): Maher puro vs ",
                           ESPECIFICACIONES[[modelo_sel]]$nombre),
         subtitle = paste0("Temporada 2024/25 | desplazamiento medio: ",
                           round(mean(abs(dif_p)), 4),
                           " | máximo: ", round(max(abs(dif_p)), 3)),
         x = "P(H) — Maher puro",
         y = paste0("P(H) — ", ESPECIFICACIONES[[modelo_sel]]$nombre),
         color = "|Δ|") +
    theme_minimal()

  print(g_comparacion)

  cat("\n=== MAGNITUD PRÁCTICA DEL AJUSTE POR ELO ===\n")
  cat("Desplazamiento de P(H) respecto de Maher puro:\n")
  cat("  medio  :", round(mean(abs(dif_p)), 4), "\n")
  cat("  mediano:", round(median(abs(dif_p)), 4), "\n")
  cat("  máximo :", round(max(abs(dif_p)), 4), "\n")
  cat(">> Significatividad estadística y relevancia práctica son cosas\n")
  cat("   distintas: un efecto detectable puede mover las probabilidades muy\n")
  cat("   poco. Este número es el que importa para el Script 07, porque define\n")
  cat("   cuántas señales de apuesta puede llegar a cambiar.\n")
}


# --- 13.3 Log-loss por especificación: validación vs test -------------------
g_metricas <- bind_rows(
  validacion |> select(modelo, log_loss) |> mutate(muestra = paste("Validación", ultima_train)),
  oos        |> select(modelo, log_loss) |> mutate(muestra = "Test 2024/25")
) |>
  ggplot(aes(x = modelo, y = log_loss, fill = muestra)) +
  geom_col(position = "dodge", alpha = 0.85) +
  coord_cartesian(ylim = c(min(c(validacion$log_loss, oos$log_loss)) * 0.98,
                           max(c(validacion$log_loss, oos$log_loss)) * 1.02)) +
  scale_fill_manual(values = c("#2980b9", "#e74c3c")) +
  labs(title    = "Log-loss por especificación",
       subtitle = "Menor es mejor | La selección se hizo sobre validación, no sobre test",
       x = NULL, y = "Log-loss", fill = NULL) +
  theme_minimal()

print(g_metricas)


# -----------------------------------------------------------------------------
# 14. GUARDAR
# -----------------------------------------------------------------------------

# Vectores de parámetros extendidos con Ipswich imputado (los que usa el 07).
# Misma regla que dentro de predecir(), aplicada al par (train completo → test).
eq_fit_final <- fit_final$equipos
nuevos_final <- setdiff(equipos_test, eq_fit_final)
ausent_final <- setdiff(eq_fit_final, equipos_test)

pesos_imp    <- ifelse(eq_fit_final %in% ausent_final, PESO_DESCENDIDOS, 1)
ataques_fin  <- fit_final$ataques
defensas_fin <- fit_final$defensas
for (eq in nuevos_final) {
  ataques_fin[eq]  <- weighted.mean(fit_final$ataques,  pesos_imp)
  defensas_fin[eq] <- weighted.mean(fit_final$defensas, pesos_imp)
}

# Consistencia: los parámetros guardados deben reproducir exactamente las
# tasas ya calculadas en resultados_fin. Si el Script 07 recalcula λ y μ a
# partir de estos vectores, tiene que obtener los mismos números.
eta_chk <- if (fit_final$nb > 0) {
  as.numeric(
    (as.matrix(datos_test_elo[, fit_final$covs, drop = FALSE]) / ESCALA_ELO) %*%
      fit_final$beta
  )
} else 0

lambda_chk <- as.numeric(exp(
  fit_final$media_global + fit_final$ventaja_local +
    ataques_fin[datos_test_elo$home_team] +
    defensas_fin[datos_test_elo$away_team] + eta_chk
))

stopifnot(max(abs(lambda_chk - resultados_fin$lambda)) < 1e-10)
cat("\nChequeo de consistencia de parámetros guardados: OK\n")

saveRDS(
  list(
    # --- salida principal para el Script 07 ---
    resultados_fin    = resultados_fin,
    modelo_final      = modelo_sel,
    nombre_modelo     = ESPECIFICACIONES[[modelo_sel]]$nombre,
    covariables       = fit_final$covs,
    escala_elo        = ESCALA_ELO,

    media_global_fin  = fit_final$media_global,
    ventaja_local_fin = fit_final$ventaja_local,
    beta_fin          = fit_final$beta,
    se_beta_fin       = fit_final$se_beta,
    ataques_fin       = ataques_fin,
    defensas_fin      = defensas_fin,

    # --- comparación de modelos ---
    fits_train_full   = fits_full,
    fits_sin_burnin   = fits_sinbi,
    fits_validacion   = fits_val,
    tabla_insample    = insample,
    tabla_validacion  = validacion,
    tabla_oos         = oos,
    tabla_sensibilidad= sensibilidad,

    # --- contrastes pareados ---
    dm_validacion     = dm_val,
    dm_test           = dm_oos,

    # --- diagnósticos ---
    calibracion       = calib_agregada,
    calibracion_test  = calib_test,
    marcadores_bajos  = marcadores,

    criterio_seleccion = paste(
      "Tres requisitos simultáneos: (a) menor log-loss",
      "sobre validación interna (fit hasta", ultima_train, "excluida,",
      "evaluación sobre", ultima_train, "); (b) contraste pareado de",
      "Diebold-Mariano contra M0 que respalde la mejora (alfa =", ALFA_DM,
      "); (c) signo de beta estable entre las muestras de estimación con y",
      "sin", primera_temp, ". Los requisitos (a) y (b) se fijaron antes de",
      "estimar las especificaciones; el (c) se incorporó durante el análisis",
      "como exigencia estándar de estabilidad paramétrica y opera solo sobre",
      "las muestras de estimación. El test 2024/25 no participó de la selección."),

    estabilidad_beta = imap_dfr(ESPECIFICACIONES, \(e, id) {
      if (length(e$covs) == 0) return(NULL)
      tibble(modelo = id, covariable = e$covs,
             beta_train_completo = fits_full[[id]]$beta,
             se_train_completo   = fits_full[[id]]$se_beta,
             beta_sin_burnin     = fits_sinbi[[id]]$beta,
             se_sin_burnin       = fits_sinbi[[id]]$se_beta,
             signo_estable = sign(fits_full[[id]]$beta) ==
                             sign(fits_sinbi[[id]]$beta))
    }),

    advertencias = c(
      paste("Sesgo de localía heredado del Script 05: verificar calibracion$brecha",
            "para el resultado H. Si es positiva, el modelo sobreestima al local",
            "en todo el test y el Script 07 debe reportar el rendimiento por",
            "lado apostado."),
      paste("Poisson independiente: revisar marcadores_bajos. El subregistro de",
            "0-0 y 1-1 es la limitación que corrige Dixon-Coles (1997)."),
      paste("El LRT/AIC/BIC in-sample está sesgado contra el Elo porque los",
            2 * (n_equipos - 1), "parámetros de equipo ya absorben la fuerza",
            "sistemática. No usarlo como criterio decisorio."),
      paste("Identificación de beta: revisar estabilidad_beta. En M1 el",
            "coeficiente se INVIERTE de signo entre muestras de estimación",
            "(positivo con", primera_temp, ", significativamente negativo sin",
            "ella). Un beta negativo implicaría que el equipo con mayor Elo",
            "convierte menos goles, lo cual carece de sentido: el coeficiente",
            "está absorbiendo ruido por colinealidad con los parámetros de",
            "equipo. Por eso la regla de selección exige estabilidad de signo."),
      paste("Rol identificante de", primera_temp, ": en esa temporada",
            "delta_elo_ini es cero para todos los equipos, de modo que",
            "delta_elo es puro componente dinámico y no colineal con Maher.",
            "Esas 380 observaciones son las que aportan la variación que",
            "identifica beta; al excluirlas el coeficiente se vuelve",
            "inestimable."),
      paste("M3 (delta_elo_ini y delta_elo_din con coeficientes libres) exhibe",
            "el signo negativo esperado bajo colinealidad con los parámetros de",
            "Maher y generaliza peor que M0. Es evidencia empírica de la",
            "advertencia formulada en el Script 05."),
      "Ipswich imputado por media ponderada (peso 3 a descendidos), igual que en el Script 04."
    )
  ),
  "resultados_modelo_final.rds"
)

cat("\n=============================================================\n")
cat("Archivo 'resultados_modelo_final.rds' guardado → usar en Script 07.\n")
cat("Modelo final:", modelo_sel, "—", ESPECIFICACIONES[[modelo_sel]]$nombre, "\n")
if (fit_final$nb > 0) {
  cat("β (por", ESCALA_ELO, "pts de ΔElo):",
      paste(fit_final$covs, "=", round(fit_final$beta, 4), collapse = " | "), "\n")
}
cat("\nRevisar antes de avanzar:\n")
cat("  - calibracion       → brecha en H (sesgo de localía)\n")
cat("  - marcadores_bajos  → subregistro de empates (Dixon-Coles)\n")
cat("  - tabla_validacion vs tabla_oos → consistencia del ranking\n")
cat("=============================================================\n")
