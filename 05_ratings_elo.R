# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 05: Construcción de ratings Elo históricos
# =============================================================================
# Depende de: Script 01 (datos_train.rds, datos_test.rds)
# Opcional:   Script 04 (resultados_maher.rds) para el diagnóstico 8.4
#
# Se construyen ratings Elo para todos los equipos procesando los partidos
# en orden cronológico. Para cada partido se guarda el Elo PRE-PARTIDO
# (antes de la actualización), que es el que luego entra como covariable
# en el modelo Maher+Elo del Script 06.
#
# -----------------------------------------------------------------------------
# FORMULACIÓN
# -----------------------------------------------------------------------------
# Resultado esperado del local (Hvattum & Arntzen, 2010):
#
#     E_H = 1 / (1 + c^((ℓ_A − ℓ_H) / d)),    con c = 10, d = 400
#
# En el código se escribe de forma algebraicamente idéntica, ordenando la
# diferencia como "local − visitante" (que es la convención del resto del
# script) y compensando con un signo menos en el exponente:
#
#     E = 1 / (1 + 10^(-(ℓ_H + h − ℓ_A) / 400))
#
# porque  -(ℓ_H − ℓ_A) = (ℓ_A − ℓ_H).  Es la misma fórmula.
#
# El término h (ventaja de localía) es una EXTENSIÓN respecto de la fórmula
# base del paper: H&A introducen el rating Elo como covariable de un modelo
# logit ordenado, y la ventaja de localía queda absorbida por los interceptos
# de esa regresión. Aquí se incorpora h dentro de la actualización para que
# el Elo sea un pronosticador autocontenido; sin h, los equipos locales
# reciben en promedio actualizaciones positivas y los visitantes negativas,
# introduciendo ruido en el rating. Se documenta como decisión propia.
#
# -----------------------------------------------------------------------------
# NOTA SOBRE AUSENCIA DE DATA LEAKAGE
# -----------------------------------------------------------------------------
# El Elo se construye sobre train + test unidos. Esto NO constituye leakage:
# el rating es estrictamente secuencial y para cada partido se registra el
# valor PREVIO a su disputa, por lo que cada observación del test usa
# únicamente información anterior a esa fecha.
# La calibración de parámetros (sección 4) se realiza exclusivamente sobre
# train. Los diagnósticos ex-post de la sección 8.5 se calculan sobre test
# pero NO retroalimentan ningún parámetro: son descriptivos.
#
# -----------------------------------------------------------------------------
# DECISIONES METODOLÓGICAS
# -----------------------------------------------------------------------------
# 1. Factor K variable:      K = K0 * (1 + δ), δ = |goles_local − goles_visit|
# 2. Ventaja de localía:     h sumada al Elo local en el resultado esperado
# 3. Regresión a la media:   entre temporadas, elo <- 1500 + (1−κ)*(elo−1500)
# 4. Equipos ascendidos:     inicializados en el Elo medio de los descendidos
#                            de la temporada anterior (post-regresión).
#                            Evita la inconsistencia de que un ascendido sin
#                            historial (Ipswich) arranque en 1500 mientras un
#                            ascendido con historial (Leicester, Southampton)
#                            arrastra un Elo congelado de dos años atrás.
# 5. PARÁMETROS FINALES — criterio de parsimonia:
#    Se corre un grid search completo sobre train con fines diagnósticos, pero
#    la superficie de error resulta prácticamente plana (ver sección 4). Dado
#    que la ganancia del óptimo del grid sobre los valores de la literatura es
#    marginal y proviene casi por completo de la regresión a la media, se
#    ADOPTAN h = 80 y K0 = 10 (Hvattum & Arntzen, 2010) y se calibra
#    únicamente κ, que no tiene referencia en el paper.
# =============================================================================

library(tidyverse)

set.seed(904158)

# -----------------------------------------------------------------------------
# 1. CARGAR DATOS Y ORDENAR CRONOLÓGICAMENTE
# -----------------------------------------------------------------------------

datos_train <- readRDS("datos_train.rds")
datos_test  <- readRDS("datos_test.rds")

datos_todos <- bind_rows(datos_train, datos_test) |>
  mutate(date = as.Date(date, format = "%d/%m/%Y")) |>
  arrange(date)

cat("Total de partidos (todas las temporadas):", nrow(datos_todos), "\n")
cat("Período:", format(min(datos_todos$date)), "→",
    format(max(datos_todos$date)), "\n")

# --- Chequeos de integridad -------------------------------------------------
stopifnot(!any(is.na(datos_todos$date)))
stopifnot(!any(is.na(datos_todos$home_team)), !any(is.na(datos_todos$away_team)))
stopifnot(!any(is.na(datos_todos$fthg)),      !any(is.na(datos_todos$ftag)))
stopifnot(all(datos_todos$ftr %in% c("H", "D", "A")))


# -----------------------------------------------------------------------------
# 2. ESTRUCTURAS AUXILIARES (índices enteros para velocidad)
# -----------------------------------------------------------------------------

equipos_todos <- sort(unique(c(datos_todos$home_team, datos_todos$away_team)))
n_equipos     <- length(equipos_todos)

temporadas    <- unique(datos_todos$temporada)   # en orden cronológico
n_temporadas  <- length(temporadas)

idx_local     <- match(datos_todos$home_team, equipos_todos)
idx_visitante <- match(datos_todos$away_team, equipos_todos)
idx_temporada <- match(datos_todos$temporada, temporadas)

stopifnot(!is.unsorted(idx_temporada))   # las temporadas no se solapan

s_observado <- dplyr::case_when(
  datos_todos$ftr == "H" ~ 1,
  datos_todos$ftr == "D" ~ 0.5,
  datos_todos$ftr == "A" ~ 0
)

delta_goles <- abs(datos_todos$fthg - datos_todos$ftag)

# Matriz lógica: ¿el equipo j juega la temporada i?
juega <- matrix(FALSE, nrow = n_temporadas, ncol = n_equipos)
for (t in seq_len(n_temporadas)) {
  eq_t <- unique(c(idx_local[idx_temporada == t], idx_visitante[idx_temporada == t]))
  juega[t, eq_t] <- TRUE
}

cat("\nEquipos distintos en la muestra:", n_equipos, "\n")
cat("Temporadas:", paste(temporadas, collapse = ", "), "\n")

cat("\n--- Movimientos entre temporadas ---\n")
for (t in 2:n_temporadas) {
  asc  <- equipos_todos[!juega[t - 1, ] &  juega[t, ]]
  desc <- equipos_todos[ juega[t - 1, ] & !juega[t, ]]
  cat(temporadas[t], "| Ascienden:", paste(asc, collapse = ", "),
      "| Descienden:", paste(desc, collapse = ", "), "\n")
}


# -----------------------------------------------------------------------------
# 3. FUNCIÓN PRINCIPAL DE CONSTRUCCIÓN DEL ELO
# -----------------------------------------------------------------------------
# Devuelve:
#   pre_l,  pre_v  : Elo ANTES del partido  (covariable del Script 06)
#   post_l, post_v : Elo DESPUÉS del partido (para el ranking histórico)
#   elo            : vector de Elo al cierre de la muestra
#   snap           : matriz (temporada x equipo) con el Elo al INICIO de cada
#                    temporada, post-regresión y post-ajuste de ascendidos

ELO_INICIAL <- 1500

construir_elo <- function(idx_l, idx_v, s_obs, dgoles, idx_temp, juega,
                          ventaja, k0, regresion, elo_inicial = ELO_INICIAL) {

  n     <- length(idx_l)
  n_eq  <- ncol(juega)
  n_tp  <- nrow(juega)

  elo    <- rep(elo_inicial, n_eq)
  pre_l  <- numeric(n); pre_v  <- numeric(n)
  post_l <- numeric(n); post_v <- numeric(n)
  snap   <- matrix(NA_real_, nrow = n_tp, ncol = n_eq)

  t_act <- idx_temp[1]
  snap[t_act, ] <- elo

  for (i in seq_len(n)) {

    # --- Corte de temporada -------------------------------------------------
    if (idx_temp[i] != t_act) {
      t_prev <- t_act
      t_act  <- idx_temp[i]

      # (a) Regresión a la media: la fuerza no persiste al 100% entre
      #     temporadas (mercado de pases, cambios de DT, plantel).
      elo <- elo_inicial + (1 - regresion) * (elo - elo_inicial)

      # (b) Ascendidos = Elo medio de los descendidos de la temporada previa.
      #     Trato idéntico para debutantes (Ipswich) y para los que vuelven
      #     tras una ausencia (Leicester, Southampton), descartando en estos
      #     últimos el Elo congelado de temporadas anteriores.
      desc <- which( juega[t_prev, ] & !juega[t_act, ])
      asc  <- which(!juega[t_prev, ] &  juega[t_act, ])
      if (length(asc) > 0 && length(desc) > 0) {
        elo[asc] <- mean(elo[desc])
      }

      snap[t_act, ] <- elo
    }

    # --- Partido ------------------------------------------------------------
    a <- idx_l[i]; b <- idx_v[i]
    ea <- elo[a];  eb <- elo[b]

    pre_l[i] <- ea
    pre_v[i] <- eb

    # E_H = 1 / (1 + 10^((ℓ_A − ℓ_H) / 400))  escrito como local − visitante
    E <- 1 / (1 + 10^(-(ea + ventaja - eb) / 400))

    K      <- k0 * (1 + dgoles[i])          # factor K variable
    ajuste <- K * (s_obs[i] - E)            # actualización suma cero

    elo[a] <- ea + ajuste
    elo[b] <- eb - ajuste

    post_l[i] <- elo[a]
    post_v[i] <- elo[b]
  }

  list(pre_l = pre_l, pre_v = pre_v,
       post_l = post_l, post_v = post_v,
       elo = elo, snap = snap)
}

# Resultado esperado, como función auxiliar reutilizable
esperado_elo <- function(el, ev, ventaja) {
  1 / (1 + 10^(-(el + ventaja - ev) / 400))
}


# -----------------------------------------------------------------------------
# 4. CALIBRACIÓN SOBRE TRAIN
# -----------------------------------------------------------------------------
# Criterio: error cuadrático medio del pronóstico Elo, mean((S − E)^2),
# equivalente al Brier score sobre la escala [0, 1].
#
# Se excluyen de la EVALUACIÓN los partidos de la primera temporada (burn-in):
# con todos los equipos en 1500 esos pronósticos son ruido. Los partidos
# igualmente se PROCESAN para construir los ratings.

es_train  <- datos_todos$temporada != "2024/25"

idx_l_tr  <- idx_local[es_train]
idx_v_tr  <- idx_visitante[es_train]
s_tr      <- s_observado[es_train]
dg_tr     <- delta_goles[es_train]
idx_tp_tr <- idx_temporada[es_train]
juega_tr  <- juega[sort(unique(idx_tp_tr)), , drop = FALSE]

eval_tr <- idx_tp_tr != min(idx_tp_tr)

cat("\nPartidos de train:", sum(es_train),
    "| usados para calibrar:", sum(eval_tr),
    "| burn-in excluido:", sum(!eval_tr), "\n")

error_elo <- function(ventaja, k0, regresion) {
  res_tr <- construir_elo(idx_l_tr, idx_v_tr, s_tr, dg_tr, idx_tp_tr, juega_tr,
                     ventaja = ventaja, k0 = k0, regresion = regresion)
  E <- esperado_elo(res_tr$pre_l, res_tr$pre_v, ventaja)
  mean((s_tr[eval_tr] - E[eval_tr])^2)
}

# --- 4.1 Grid search completo (diagnóstico) ---------------------------------

grilla <- expand_grid(
  ventaja   = seq(0, 140, by = 10),
  k0        = seq(4, 40, by = 2),
  regresion = c(0, 0.10, 0.20, 0.30, 0.40)
)

cat("\nCalibrando", nrow(grilla), "combinaciones sobre train...\n")
t0 <- Sys.time()
grilla$error <- pmap_dbl(grilla, error_elo)
cat("Completado en",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "segundos\n")

mejor_grid <- grilla |> slice_min(error, n = 1)

cat("\n--- Top 10 combinaciones del grid ---\n")
print(grilla |> slice_min(error, n = 10))


# --- 4.2 Parámetros adoptados: literatura + κ calibrado ---------------------
# H&A (2010) no contempla regresión a la media entre temporadas, así que κ es
# el único parámetro sin referencia bibliográfica y el único que se estima.

VENTAJA_LOCAL <- 80   # Hvattum & Arntzen (2010)
K0            <- 10   # Hvattum & Arntzen (2010)

grid_r <- tibble(regresion = seq(0, 0.50, by = 0.05)) |>
  mutate(error = map_dbl(regresion, \(k) error_elo(VENTAJA_LOCAL, K0, k)))

cat("\n--- Calibración de la regresión a la media (con h=80, K0=10) ---\n")
print(grid_r |> mutate(error = round(error, 5)))

REGRESION <- grid_r |> slice_min(error, n = 1) |> pull(regresion)


# --- 4.3 Descomposición de la ganancia --------------------------------------
# Aísla cuánto aporta cada decisión respecto de la especificación base.

err_base    <- error_elo(80, 10, 0)                        # H&A puro
err_adoptado<- error_elo(VENTAJA_LOCAL, K0, REGRESION)     # H&A + regresión
err_optimo  <- mejor_grid$error                            # óptimo del grid

cat("\n--- Descomposición de la ganancia (ECM sobre train) ---\n")
cat("(a) H&A puro (h=80, K0=10, κ=0)          :", round(err_base, 5), "\n")
cat("(b) H&A + regresión (κ =", REGRESION, ")           :", round(err_adoptado, 5), "\n")
cat("(c) Óptimo del grid (h =", mejor_grid$ventaja,
    ", K0 =", mejor_grid$k0, ", κ =", mejor_grid$regresion, "):",
    round(err_optimo, 5), "\n")
cat("\nGanancia por regresión a la media  (a→b):",
    round(err_base - err_adoptado, 6),
    "(", round(100 * (err_base - err_adoptado) / err_base, 2), "% )\n")
cat("Ganancia adicional por reoptimizar h y K0 (b→c):",
    round(err_adoptado - err_optimo, 6),
    "(", round(100 * (err_adoptado - err_optimo) / err_base, 2), "% )\n")

cat("\n>> DECISIÓN: se adoptan h =", VENTAJA_LOCAL, "y K0 =", K0,
    "(literatura), κ =", REGRESION, "(calibrado).\n")
cat(">> La reoptimización de h y K0 aporta una fracción marginal del error y\n")
cat("   la superficie es plana (ver top 10 y gráfico 10.5): se prioriza\n")
cat("   parsimonia y comparabilidad con la literatura.\n")


# -----------------------------------------------------------------------------
# 5. CONSTRUCCIÓN FINAL DEL ELO (TRAIN + TEST)
# -----------------------------------------------------------------------------

cat("\nProcesando", nrow(datos_todos), "partidos con los parámetros adoptados...\n")

res <- construir_elo(idx_local, idx_visitante, s_observado, delta_goles,
                     idx_temporada, juega,
                     ventaja   = VENTAJA_LOCAL,
                     k0        = K0,
                     regresion = REGRESION)

stopifnot(!any(is.na(res$pre_l)), !any(is.na(res$pre_v)))
cat("Procesamiento completado ✓\n")

elo <- setNames(res$elo, equipos_todos)


# -----------------------------------------------------------------------------
# 6. AGREGAR ELO AL DATASET
# -----------------------------------------------------------------------------
# Descomposición del ΔElo en dos componentes causales:
#
#   delta_elo_ini = diferencial al INICIO de la temporada
#                   → componente ESTÁTICO. Fuerza relativa que los parámetros
#                     de ataque/defensa de Maher (Script 04) ya capturan por
#                     construcción, al ser fijos dentro de la temporada.
#
#   delta_elo_din = delta_elo − delta_elo_ini
#                   → componente DINÁMICO. Movimiento acumulado desde el
#                     arranque de la temporada; información que Maher NO puede
#                     capturar. Es el aporte específico del Elo.
#
# Ambos son causales: el snapshot inicial se fija antes del primer partido de
# la temporada, sin usar información futura.

elo_local_ini     <- res$snap[cbind(idx_temporada, idx_local)]
elo_visitante_ini <- res$snap[cbind(idx_temporada, idx_visitante)]

datos_todos <- datos_todos |>
  mutate(
    elo_local         = res$pre_l,
    elo_visitante     = res$pre_v,
    delta_elo         = elo_local - elo_visitante,

    elo_local_post    = res$post_l,
    elo_visitante_post= res$post_v,

    elo_local_ini     = elo_local_ini,
    elo_visitante_ini = elo_visitante_ini,
    delta_elo_ini     = elo_local_ini - elo_visitante_ini,
    delta_elo_din     = delta_elo - delta_elo_ini,

    elo_esperado      = esperado_elo(elo_local, elo_visitante, VENTAJA_LOCAL),
    s_observado       = s_observado
  )


# -----------------------------------------------------------------------------
# 7. SEPARAR EN TRAIN Y TEST
# -----------------------------------------------------------------------------

datos_train_elo <- datos_todos |> filter(temporada != "2024/25")
datos_test_elo  <- datos_todos |> filter(temporada == "2024/25")

cat("\nPartidos train con Elo:", nrow(datos_train_elo), "\n")
cat("Partidos test  con Elo:", nrow(datos_test_elo),  "\n")


# -----------------------------------------------------------------------------
# 8. DIAGNÓSTICOS
# -----------------------------------------------------------------------------

primera_temp <- temporadas[1]

# --- 8.1 Burn-in: dispersión del ΔElo por temporada -------------------------

cat("\n--- 8.1 Burn-in: dispersión del ΔElo por temporada ---\n")
datos_todos |>
  group_by(temporada) |>
  summarise(
    n            = n(),
    sd_delta     = round(sd(delta_elo), 1),
    sd_delta_din = round(sd(delta_elo_din), 1),
    rango        = paste0(round(min(delta_elo)), " a ", round(max(delta_elo)))
  ) |>
  print()

cat("\n>> En", primera_temp, "el ΔElo inicial es idénticamente 0 (todos en 1500):\n")
cat("   esa temporada aporta señal atenuada. Evaluar en el Script 06 estimar\n")
cat("   con y sin ella y comparar.\n")


# --- 8.2 Descomposición estático / dinámico ---------------------------------
# Se excluye la primera temporada: sd(delta_elo_ini) = 0 y cor() no está
# definida (todos los equipos arrancan en el mismo valor).

cat("\n--- 8.2 Descomposición del ΔElo (estático vs dinámico) ---\n")

diag_desc <- datos_todos |>
  filter(temporada != primera_temp) |>
  group_by(temporada) |>
  summarise(
    sd_total    = sd(delta_elo),
    sd_estatico = sd(delta_elo_ini),
    sd_dinamico = sd(delta_elo_din),
    cor_tot_est = cor(delta_elo, delta_elo_ini),
    r2_estatico = cor(delta_elo, delta_elo_ini)^2,
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

print(diag_desc)

cor_test <- cor(datos_test_elo$delta_elo, datos_test_elo$delta_elo_ini)
cat("\nEn test (24/25): cor(ΔElo, ΔElo inicial) =", round(cor_test, 3),
    "→ R² =", round(cor_test^2, 3), "\n")
cat("Varianza del ΔElo no explicada por el nivel inicial:",
    round(100 * (1 - cor_test^2), 1), "%\n")
cat("Desvío del componente dinámico:",
    round(sd(datos_test_elo$delta_elo_din), 1), "puntos Elo",
    "≈", round(100 * (esperado_elo(sd(datos_test_elo$delta_elo_din), 0, 0) - 0.5), 1),
    "puntos porcentuales de probabilidad\n")
cat(">> Un R² alto NO implica que el componente dinámico sea irrelevante:\n")
cat("   indica que el ΔElo total está dominado por el nivel de arranque.\n")
cat(">> Especificación recomendada para el Script 06: delta_elo_din.\n")


# --- 8.3 Calibración del pronóstico Elo -------------------------------------

cat("\n--- 8.3 Calibración del pronóstico Elo en test (24/25) ---\n")
datos_test_elo |>
  mutate(bin = cut(elo_esperado, breaks = seq(0, 1, by = 0.1),
                   include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(n = n(),
            esperado  = round(mean(elo_esperado), 3),
            observado = round(mean(s_observado),  3),
            .groups = "drop") |>
  print()

cat("\nECM del pronóstico Elo — train (sin burn-in):",
    round(mean((datos_train_elo$s_observado[datos_train_elo$temporada != primera_temp] -
                datos_train_elo$elo_esperado[datos_train_elo$temporada != primera_temp])^2), 5), "\n")
cat("ECM del pronóstico Elo — test (24/25)       :",
    round(mean((datos_test_elo$s_observado - datos_test_elo$elo_esperado)^2), 5), "\n")


# --- 8.4 Colinealidad con los parámetros de Maher (Script 04) ---------------
# ACLARACIÓN: la colinealidad NO invalida la comparación de modelos anidados
# del Script 06 (LRT / AIC / log-loss fuera de muestra). Afecta la lectura del
# coeficiente individual del ΔElo y la estabilidad numérica del optimizador.
# Sirve para interpretar el resultado, no para decidir si el Script 06 vale.

ARCHIVO_MAHER <- "resultados_maher.rds"

cat("\n--- 8.4 Colinealidad con parámetros de Maher ---\n")

if (file.exists(ARCHIVO_MAHER)) {

  maher <- readRDS(ARCHIVO_MAHER)
  at <- maher$ataques
  df <- maher$defensas

  faltantes <- setdiff(equipos_todos, names(at))
  if (length(faltantes) > 0) {
    cat("AVISO: equipos sin parámetro de Maher:", paste(faltantes, collapse = ", "), "\n")
  }

  # Fuerza relativa implícita en Maher:
  #   ataque propio − ataque rival  +  debilidad defensiva rival − propia
  datos_todos <- datos_todos |>
    mutate(
      fuerza_maher = as.numeric(
        (at[home_team] - at[away_team]) + (df[away_team] - df[home_team])
      )
    )

  datos_todos |>
    filter(temporada != primera_temp, !is.na(fuerza_maher)) |>
    group_by(temporada) |>
    summarise(
      cor_total    = round(cor(delta_elo,     fuerza_maher), 3),
      cor_estatico = round(cor(delta_elo_ini, fuerza_maher), 3),
      cor_dinamico = round(cor(delta_elo_din, fuerza_maher), 3),
      .groups = "drop"
    ) |>
    print()

  cat("\n>> ATENCIÓN AL LEER ESTA TABLA:\n")
  cat("   Los parámetros de Maher se estimaron sobre TODO el train, de modo que\n")
  cat("   para 2022/23 y 2023/24 el modelo ya vio esos partidos y absorbió parte\n")
  cat("   del movimiento intra-temporada. Las cor_dinamico de esas filas están\n")
  cat("   contaminadas y aparecen infladas.\n")
  cat("   La ÚNICA medición limpia es la de 2024/25: Maher no vio esa temporada.\n")
  cat(">> Si cor_dinamico en 2024/25 es baja, delta_elo_din es aproximadamente\n")
  cat("   ortogonal a Maher fuera de muestra y es la covariable correcta para\n")
  cat("   el Script 06.\n")

  # Se resincroniza el split tras agregar la columna
  datos_train_elo <- datos_todos |> filter(temporada != "2024/25")
  datos_test_elo  <- datos_todos |> filter(temporada == "2024/25")

} else {
  cat("No se encontró '", ARCHIVO_MAHER, "'. Diagnóstico omitido.\n", sep = "")
}


# --- 8.5 DERIVA DE LA VENTAJA DE LOCALÍA ------------------------------------
# HALLAZGO: el pronóstico Elo sobreestima al local en la 24/25.
#
# ADVERTENCIA SOBRE LA FUERZA DE LA EVIDENCIA:
# El Script 04 detecta lo mismo por otra vía (brecha_ventaja_pct), pero eso NO
# constituye confirmación independiente: ambos modelos se calculan sobre LOS
# MISMOS 380 partidos. Si la 24/25 tuvo pocas victorias locales por azar
# muestral, los dos modelos lo reflejarían igual. Son dos mediciones del mismo
# dato, no dos fuentes de evidencia. Debe presentarse como evidencia
# sugestiva, no concluyente.
#
# TEST CORRECTO: no alcanza con contrastar el sesgo del test contra cero,
# porque h = 80 ya produce sesgo en train. Hay que comparar el sesgo del TEST
# contra el sesgo del TRAIN (diferencia en diferencias).
#
# CONSECUENCIA PARA LOS SCRIPTS DE APUESTAS: ambos modelos van a asignar
# probabilidad excesiva al local durante el test, generando señales de "valor"
# potencialmente espurias en apuestas al local. Controlar reportando el
# rendimiento por separado según el lado apostado.
#
# Estos cálculos son DESCRIPTIVOS: no retroalimentan ningún parámetro.

cat("\n--- 8.5 Deriva de la ventaja de localía ---\n")

# Ventaja de localía implícita: valor de h que iguala E medio con S medio
ventaja_implicita <- function(el, ev, s) {
  f <- function(v) mean(esperado_elo(el, ev, v)) - mean(s)
  uniroot(f, interval = c(-400, 400))$root
}

deriva <- datos_todos |>
  group_by(temporada) |>
  summarise(
    n            = n(),
    n_local      = sum(ftr == "H"),
    pct_local    = round(100 * mean(ftr == "H"), 1),
    pct_empate   = round(100 * mean(ftr == "D"), 1),
    pct_visita   = round(100 * mean(ftr == "A"), 1),
    s_observado  = round(mean(s_observado), 3),
    e_esperado   = round(mean(elo_esperado), 3),
    brecha       = round(mean(elo_esperado) - mean(s_observado), 3),
    h_implicita  = round(ventaja_implicita(elo_local, elo_visitante, s_observado), 1),
    .groups = "drop"
  ) |>
  # IC binomial del % de victorias locales: contextualiza la variabilidad
  mutate(ic95_pct_local = map2_chr(n_local, n, \(x, nn) {
    ci <- binom.test(x, nn)$conf.int
    paste0("[", round(100 * ci[1], 1), ", ", round(100 * ci[2], 1), "]")
  }))

print(deriva, width = Inf)

cat("\nVentaja de localía asumida: h =", VENTAJA_LOCAL, "puntos Elo\n")
cat("Ventaja implícita por temporada:",
    paste(temporadas, "=", deriva$h_implicita, collapse = " | "), "\n")
cat(">> La serie de h implícita es volátil: los IC del % de victorias locales\n")
cat("   se solapan entre temporadas. No sobreinterpretar un único valor.\n")

# --- Test correcto: sesgo del TEST contra sesgo del TRAIN -------------------
# El sesgo del test contra cero NO es informativo, porque h = 80 ya produce
# sesgo en train. Lo que interesa es si el sesgo CAMBIÓ.

sesgo_train <- datos_train_elo |>
  filter(temporada != primera_temp) |>
  mutate(sesgo = s_observado - elo_esperado) |>
  pull(sesgo)

sesgo_test <- datos_test_elo |>
  mutate(sesgo = s_observado - elo_esperado) |>
  pull(sesgo)

cat("\n--- Sesgo medio (S − E) ---\n")
cat("Train sin burn-in (", length(sesgo_train), "partidos):",
    round(mean(sesgo_train), 4), "\n")
cat("Test 24/25        (", length(sesgo_test),  "partidos):",
    round(mean(sesgo_test),  4), "\n")

tt_ref <- t.test(sesgo_test)                    # referencia: sesgo vs cero
tt_dif <- t.test(sesgo_test, sesgo_train)       # test relevante

cat("\n(a) Referencia — sesgo del test contra cero:\n")
cat("    estimación:", round(unname(tt_ref$estimate), 4),
    "| p =", format.pval(tt_ref$p.value, digits = 3),
    "  <-- NO es el test de interés\n")

dif <- unname(tt_dif$estimate[1] - tt_dif$estimate[2])
cat("\n(b) TEST RELEVANTE — sesgo del test contra el sesgo del train:\n")
cat("    diferencia:", round(dif, 4),
    "| IC 95%: [", round(tt_dif$conf.int[1], 4), ",",
    round(tt_dif$conf.int[2], 4), "]",
    "| p =", format.pval(tt_dif$p.value, digits = 3), "\n")

if (tt_dif$p.value < 0.05) {
  cat("\n>> El sesgo del test difiere del de train al 5%. Evidencia de que la\n")
  cat("   ventaja de localía cayó en la 24/25 — pero verificar cuán cerca del\n")
  cat("   umbral está el p-valor antes de afirmarlo con firmeza.\n")
} else {
  cat("\n>> No se rechaza que el sesgo del test sea igual al de train: la\n")
  cat("   sobreestimación del local es atribuible a h = 80, no a un cambio\n")
  cat("   estructural en la 24/25.\n")
}

# --- Robustez: ¿el resultado depende de excluir el burn-in? -----------------
# La línea de base del test usa train SIN 2021/22. Esa temporada tuvo una
# ventaja de localía baja (ver h_implicita), de modo que incluirla desplaza
# la referencia y puede revertir la conclusión. El jurado puede preguntarlo.

sesgo_train_full <- datos_train_elo |>
  mutate(sesgo = s_observado - elo_esperado) |>
  pull(sesgo)

tt_dif_full <- t.test(sesgo_test, sesgo_train_full)
dif_full <- unname(tt_dif_full$estimate[1] - tt_dif_full$estimate[2])

cat("\n--- Robustez del test (b) ---\n")
cat("Línea de base SIN 2021/22 (", length(sesgo_train), "part.): dif =",
    round(dif, 4), "| p =", format.pval(tt_dif$p.value, digits = 3), "\n")
cat("Línea de base CON 2021/22 (", length(sesgo_train_full), "part.): dif =",
    round(dif_full, 4), "| p =", format.pval(tt_dif_full$p.value, digits = 3), "\n")

if (xor(tt_dif$p.value < 0.05, tt_dif_full$p.value < 0.05)) {
  cat("\n>> LA CONCLUSIÓN CAMBIA SEGÚN SE INCLUYA O NO EL BURN-IN.\n")
  cat("   La exclusión de 2021/22 está justificada metodológicamente (ratings\n")
  cat("   inmaduros), pero el resultado es frágil. Reportar AMBAS versiones y\n")
  cat("   presentar la deriva de localía como indicio, no como hallazgo firme.\n")
} else {
  cat("\n>> La conclusión es estable ante la inclusión del burn-in.\n")
}

# --- Lectura recomendada de la serie de h implícita -------------------------
cat("\n>> LECTURA PARA LA REDACCIÓN:\n")
cat("   La serie de h implícita (", paste(deriva$h_implicita, collapse = ", "),
    ") no describe\n", sep = "")
cat("   una tendencia decreciente sino VOLATILIDAD entre temporadas: la 2021/22\n")
cat("   es casi tan baja como la 2024/25. La afirmación defendible es que\n")
cat("   h = 80 —calibrado por H&A con datos 1993-2008— resulta alto para el\n")
cat("   período analizado, no que exista un quiebre estructural en la 24/25.\n")

# Cotejo con el Script 04
if (exists("maher")) {
  cat("\n--- Cotejo con el Script 04 (Maher) ---\n")
  cat("vent_modelo:", maher$vent_modelo,
      "| vent_test:", maher$vent_test,
      "| brecha:", maher$brecha_ventaja_pct, "%\n")
  if (!is.null(maher$ventaja_por_temporada)) {
    print(maher$ventaja_por_temporada, width = Inf)
  }
  cat(">> Maher y Elo apuntan en la misma dirección, PERO se calculan sobre los\n")
  cat("   mismos 380 partidos: es la misma evidencia medida dos veces, no dos\n")
  cat("   confirmaciones independientes.\n")
}


# -----------------------------------------------------------------------------
# 9. RANKING ELO
# -----------------------------------------------------------------------------
# El Elo de un equipo ausente de la Premier sigue siendo alcanzado por la
# regresión a la media en cada corte de temporada, de modo que el valor final
# del vector 'elo' NO coincide con su rating tras el último partido disputado.
# Para el ranking histórico se usa el Elo POST-PARTIDO de su último encuentro.

elo_ultimo <- datos_todos |>
  select(date, temporada, home_team, away_team, elo_local_post, elo_visitante_post) |>
  pivot_longer(
    cols = c(home_team, away_team),
    names_to = "condicion", values_to = "equipo"
  ) |>
  mutate(elo_post = if_else(condicion == "home_team",
                            elo_local_post, elo_visitante_post)) |>
  arrange(equipo, date) |>
  group_by(equipo) |>
  slice_tail(n = 1) |>
  ungroup() |>
  select(equipo,
         elo_ultimo_partido = elo_post,
         ultima_fecha       = date,
         ultima_temporada   = temporada)

elo_final <- tibble(
  equipo    = equipos_todos,
  elo_cierre = round(as.numeric(elo), 1)   # con regresiones posteriores aplicadas
) |>
  left_join(elo_ultimo, by = "equipo") |>
  mutate(elo_ultimo_partido = round(elo_ultimo_partido, 1)) |>
  arrange(desc(elo_ultimo_partido))

equipos_test <- sort(unique(c(datos_test$home_team, datos_test$away_team)))

cat("\n--- Elo tras el último partido disputado por cada equipo ---\n")
cat("(ultima_temporada indica a qué momento corresponde el rating;\n")
cat(" elo_cierre incluye las regresiones a la media posteriores y solo es\n")
cat(" comparable entre equipos que llegaron a la 24/25)\n")
print(elo_final, n = n_equipos)

cat("\n--- Elo al cierre de la 24/25 (solo equipos que la disputaron) ---\n")
elo_final_2425 <- elo_final |>
  filter(equipo %in% equipos_test) |>
  select(equipo, elo = elo_ultimo_partido)
print(elo_final_2425, n = nrow(elo_final_2425))

cat("\n--- Distribución del ΔElo en test (24/25) ---\n")
datos_test_elo |>
  summarise(
    media  = round(mean(delta_elo), 1),
    desvio = round(sd(delta_elo),   1),
    min    = round(min(delta_elo),  1),
    max    = round(max(delta_elo),  1)
  ) |>
  print()


# -----------------------------------------------------------------------------
# 10. GRÁFICOS
# -----------------------------------------------------------------------------

subtitulo_params <- paste0("h = ", VENTAJA_LOCAL, " | K0 = ", K0,
                           " | regresión = ", REGRESION)

# --- 10.1 Evolución del Elo -------------------------------------------------
# Series agrupadas por (equipo, temporada): la línea no conecta temporadas
# distintas, de modo que se ven los saltos por regresión a la media y los
# huecos de los equipos ausentes.

equipos_seleccion <- c("Man City", "Man United", "Liverpool",
                       "Chelsea", "Arsenal", "Tottenham", "Ipswich", "Burnley")

g_evolucion <- datos_todos |>
  pivot_longer(cols = c(elo_local, elo_visitante),
               names_to = "condicion", values_to = "elo_valor") |>
  mutate(equipo = if_else(condicion == "elo_local", home_team, away_team)) |>
  filter(equipo %in% equipos_seleccion) |>
  ggplot(aes(x = date, y = elo_valor, color = equipo,
             group = interaction(equipo, temporada))) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 1500, linetype = "dashed", color = "gray50") +
  labs(
    title    = "Evolución del Elo — Big 6 + Ipswich + Burnley",
    subtitle = paste0("Premier League 2021/22 - 2024/25 | ", subtitulo_params,
                      "\nLíneas cortadas entre temporadas: regresión a la media y ausencias"),
    x = NULL, y = "Elo pre-partido", color = "Equipo"
  ) +
  theme_minimal()

print(g_evolucion)


# --- 10.2 Ranking Elo al cierre de la 24/25 ---------------------------------
g_ranking <- elo_final_2425 |>
  ggplot(aes(x = reorder(equipo, elo), y = elo, fill = elo > 1500)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 1500, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2980b9"),
                    guide = "none") +
  coord_flip() +
  labs(title    = "Ranking Elo al cierre de la temporada 24/25",
       subtitle = paste0("Solo equipos de la 24/25 | ", subtitulo_params),
       x = NULL, y = "Elo") +
  theme_minimal()

print(g_ranking)


# --- 10.3 Variación del Elo durante la 24/25 --------------------------------
t_2425 <- match("2024/25", temporadas)

elo_variacion <- tibble(
  equipo     = equipos_todos,
  elo_inicio = res$snap[t_2425, ]
) |>
  inner_join(elo_final_2425 |> rename(elo_fin = elo), by = "equipo") |>
  mutate(variacion = elo_fin - elo_inicio) |>
  arrange(desc(variacion))

cat("\n--- Variación del Elo en la 24/25 ---\n")
print(elo_variacion |> mutate(across(where(is.numeric), \(x) round(x, 1))), n = 20)

g_variacion <- elo_variacion |>
  ggplot(aes(x = reorder(equipo, variacion), y = variacion, fill = variacion > 0)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2980b9"),
                    guide = "none") +
  coord_flip() +
  labs(title    = "Variación del Elo durante la temporada 24/25",
       subtitle = "Desde el rating de arranque de temporada | Rojo = ganó Elo | Azul = perdió Elo",
       x = NULL, y = "Variación Elo (fin − inicio)") +
  theme_minimal()

print(g_variacion)


# --- 10.4 ΔElo según resultado observado ------------------------------------
g_boxplot <- datos_test_elo |>
  mutate(resultado = case_when(
    ftr == "H" ~ "Victoria local",
    ftr == "D" ~ "Empate",
    ftr == "A" ~ "Victoria visitante"
  )) |>
  ggplot(aes(x = resultado, y = delta_elo, fill = resultado)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Victoria local"     = "#e74c3c",
                               "Empate"             = "#95a5a6",
                               "Victoria visitante" = "#2980b9"), guide = "none") +
  labs(title    = "Diferencial de Elo según resultado — Temporada 24/25",
       subtitle = "ΔElo = Elo local − Elo visitante | Positivo = local favorito",
       x = NULL, y = "ΔElo pre-partido") +
  theme_minimal()

print(g_boxplot)


# --- 10.5 Superficie de error del grid --------------------------------------
# Documenta visualmente la planicie que justifica la elección por parsimonia.

g_grid <- grilla |>
  filter(regresion == REGRESION) |>
  ggplot(aes(x = ventaja, y = k0, fill = error)) +
  geom_tile() +
  geom_point(aes(x = VENTAJA_LOCAL, y = K0), color = "white",
             shape = 4, size = 4, stroke = 1.5, inherit.aes = FALSE) +
  geom_point(data = mejor_grid, aes(x = ventaja, y = k0), color = "cyan",
             shape = 1, size = 5, stroke = 1.2, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(title    = "Error cuadrático medio del pronóstico Elo (train)",
       subtitle = paste0("Corte en regresión = ", REGRESION,
                         " | X blanca = adoptado (H&A) | círculo = óptimo del grid"),
       x = "Ventaja de localía (puntos Elo)", y = "K0", fill = "ECM") +
  theme_minimal()

print(g_grid)


# --- 10.6 Componente dinámico vs estático -----------------------------------
g_dinamico <- datos_todos |>
  filter(temporada != primera_temp) |>
  ggplot(aes(x = delta_elo_ini, y = delta_elo_din)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ temporada) +
  labs(title    = "Componente dinámico vs estático del ΔElo",
       subtitle = "Eje Y = movimiento acumulado dentro de la temporada (lo que Maher no captura)",
       x = "ΔElo al inicio de la temporada", y = "ΔElo dinámico") +
  theme_minimal()

print(g_dinamico)


# --- 10.7 Calibración: esperado vs observado --------------------------------
# Visualiza el hallazgo de 8.5: los puntos por debajo de la diagonal indican
# sobreestimación del local.

g_calib <- datos_todos |>
  mutate(muestra = if_else(temporada == "2024/25", "Test (24/25)", "Train"),
         bin = cut(elo_esperado, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) |>
  group_by(muestra, bin) |>
  summarise(n = n(),
            esperado  = mean(elo_esperado),
            observado = mean(s_observado), .groups = "drop") |>
  ggplot(aes(x = esperado, y = observado, color = muestra, size = n)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Train" = "#2980b9", "Test (24/25)" = "#e74c3c")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title    = "Calibración del pronóstico Elo",
       subtitle = "Puntos bajo la diagonal = el modelo sobreestima al local",
       x = "Probabilidad esperada (Elo)", y = "Frecuencia observada",
       color = NULL, size = "n partidos") +
  theme_minimal()

print(g_calib)


# -----------------------------------------------------------------------------
# 11. GUARDAR
# -----------------------------------------------------------------------------

parametros_elo <- list(
  elo_inicial      = ELO_INICIAL,
  ventaja_local    = VENTAJA_LOCAL,
  k0               = K0,
  regresion_media  = REGRESION,
  fuente_h_k0      = "Hvattum & Arntzen (2010)",
  fuente_regresion = "Calibrado sobre train (ECM), sin burn-in",
  ecm_base_HA      = err_base,
  ecm_adoptado     = err_adoptado,
  ecm_optimo_grid  = err_optimo,
  optimo_grid      = mejor_grid,
  criterio         = "ECM del pronóstico Elo sobre train, excluyendo 2021/22 (burn-in)",
  nota_parsimonia  = paste("Se adoptan los valores de la literatura para h y K0:",
                           "la superficie de error es plana y la ganancia por",
                           "reoptimizarlos es marginal frente al riesgo de sobreajuste.")
)

diagnostico_localia <- list(
  tabla_deriva     = deriva,
  sesgo_train      = mean(sesgo_train),
  sesgo_test       = mean(sesgo_test),
  test_vs_cero     = list(estimate = unname(tt_ref$estimate),
                          p_value  = tt_ref$p.value,
                          nota     = "No es el test de interés"),
  test_vs_train    = list(diferencia = dif,
                          ic95       = as.numeric(tt_dif$conf.int),
                          p_value    = tt_dif$p.value,
                          nota       = "Test relevante: cambio del sesgo"),
  robustez_burnin  = list(diferencia = dif_full,
                          p_value    = tt_dif_full$p.value,
                          nota       = "Mismo test con 2021/22 en la línea de base"),
  interpretacion = paste("h = 80 (H&A, datos 1993-2008) resulta alto para el",
                         "período analizado. La serie de h implícita es VOLÁTIL",
                         "entre temporadas (2021/22 es casi tan baja como",
                         "2024/25): no hay tendencia decreciente, no afirmar",
                         "quiebre estructural. Evidencia SUGESTIVA, no",
                         "concluyente: Maher y Elo miden los mismos 380 partidos,",
                         "los IC del % de victorias locales se solapan y el test",
                         "es sensible a excluir el burn-in. Controlar en los",
                         "scripts de apuestas reportando resultados por lado",
                         "apostado.")
)

saveRDS(
  list(
    datos_train_elo     = datos_train_elo,
    datos_test_elo      = datos_test_elo,
    datos_todos         = datos_todos,
    elo_final           = elo_final,
    elo_final_2425      = elo_final_2425,
    elo_variacion       = elo_variacion,
    parametros_elo      = parametros_elo,
    diagnostico_localia = diagnostico_localia,
    diag_descomposicion = diag_desc,
    grilla_calib        = grilla,
    grid_regresion      = grid_r,
    snapshots           = res$snap,
    equipos             = equipos_todos,
    temporadas          = temporadas
  ),
  "resultados_elo.rds"
)

cat("\n=============================================================\n")
cat("Archivo 'resultados_elo.rds' guardado → usar en Script 06.\n")
cat("Parámetros: h =", VENTAJA_LOCAL, "(H&A) | K0 =", K0, "(H&A) | κ =",
    REGRESION, "(calibrado)\n")
cat("\nCovariables disponibles para el Script 06:\n")
cat("  delta_elo      — diferencial total\n")
cat("  delta_elo_ini  — componente estático (solapado con Maher)\n")
cat("  delta_elo_din  — componente dinámico (recomendada)\n")
cat("\nADVERTENCIA A ARRASTRAR: indicios de caída de la ventaja de localía en\n")
cat("la 24/25 (evidencia sugestiva, no concluyente). Ambos modelos tienden a\n")
cat("sobreestimar al local en test → posibles señales de valor espurias en\n")
cat("apuestas al local. Ver diagnostico_localia en el .rds.\n")
cat("=============================================================\n")
