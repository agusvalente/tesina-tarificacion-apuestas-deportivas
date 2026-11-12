# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 02: Modelo Poisson bivariado con corrección τ (parámetro ρ)
# =============================================================================
# Depende de: Script 01 (datos_train.rds, datos_test.rds)
# Su resultado se usa en el Script 03 para testear si ρ es significativo.
# =============================================================================
# install.packages("tidyverse")
# =============================================================================
#
# DECISIONES DE DISEÑO
#
#   (a) ALCANCE DEL MODELO: qué se implementa y qué no.
#       Dixon & Coles (1997) tienen DOS componentes:
#         1. la corrección τ para marcadores bajos, gobernada por ρ;
#         2. una ponderación exponencial decreciente φ(t) sobre la
#            verosimilitud, que da más peso a los partidos recientes.
#       Este script implementa SOLO la primera. Lo que se estima es, con
#       precisión, el modelo de Maher (1982) con corrección τ.
#
#       Consecuencia para la interpretación: el resultado del Script 03 NO
#       permite afirmar "Dixon-Coles no aporta sobre Maher". La afirmación
#       defendible es más acotada: la corrección τ no aporta ajuste
#       significativo sobre la especificación Poisson independiente EN ESTA
#       MUESTRA. La ausencia de ponderación temporal se declara como
#       limitación aparte.
#
#   (b) IDENTIFICACIÓN: fijar α₁ = 0 no alcanza. La parametrización tiene DOS
#       redundancias y esa restricción remueve solo una. Ver el detalle en la
#       Sección 2 y la comprobación numérica en la Sección 4.
#       Solución: normalización posterior a la estimación con DOS restricciones
#       de escala, media(ataques) = 0 y media(defensas) = 0. No se toca el
#       optimizador, así que la verosimilitud, ρ y todas las probabilidades del
#       modelo no se ven afectadas; lo único que cambia es cómo se reparte un
#       nivel constante entre media_global, ataques y defensas, y con ello la
#       interpretabilidad de los tres.
#
#   (c) Se reporta el recuento de partidos por equipo en el entrenamiento:
#       los 25 equipos NO tienen la misma cantidad de partidos y por lo tanto
#       la precisión de sus parámetros es heterogénea.
#
#   (d) La penalización del dominio inválido de τ crece con la distancia al
#       dominio en lugar de ser una constante: una constante daría gradiente
#       numérico nulo y el optimizador podría quedar clavado.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. CARGAR DATOS
# -----------------------------------------------------------------------------

datos_train <- readRDS("datos_train.rds")   # 3 temporadas → para ENTRENAR
datos_test  <- readRDS("datos_test.rds")    # 24/25        → para VALIDAR

equipos   <- sort(unique(c(datos_train$home_team, datos_train$away_team)))
n_equipos <- length(equipos)

cat("Partidos de entrenamiento:", nrow(datos_train), "\n")
cat("Equipos en entrenamiento: ", n_equipos, "\n")

# --- Partidos por equipo -----------------------------------------------------
# Sobre 3 temporadas de 20 equipos hay ascensos y descensos: aparecen 25
# equipos distintos. Un equipo presente en las tres temporadas aporta 114
# partidos para estimar sus dos parámetros; uno presente en una sola aporta 38.
# El modelo no pondera esa diferencia: los parámetros de los equipos con menos
# presencia son sensiblemente menos precisos.
# (Este recuento va al marco metodológico, junto con la cobertura de cuotas.)

partidos_por_equipo <- datos_train |>
  pivot_longer(c(home_team, away_team), values_to = "equipo") |>
  count(equipo, name = "partidos") |>
  arrange(partidos, equipo)

cat("\n--- Partidos por equipo en el entrenamiento ---\n")
print(partidos_por_equipo, n = Inf)

cat(sprintf(
  "\nRango: %d a %d partidos. Equipos con menos de 114 partidos (no jugaron\nlas tres temporadas): %d de %d.\n",
  min(partidos_por_equipo$partidos),
  max(partidos_por_equipo$partidos),
  sum(partidos_por_equipo$partidos < 114),
  n_equipos
))


# -----------------------------------------------------------------------------
# 2. FUNCIÓN DE LOG-VEROSIMILITUD CON CORRECCIÓN τ
# -----------------------------------------------------------------------------
# El modelo asume que los goles de cada equipo siguen distribuciones Poisson:
#
#   Goles_local     ~ Poisson(lambda)
#   Goles_visitante ~ Poisson(mu)
#
# Donde:
#   log(lambda) = media_global + ataque_local + defensa_visitante + ventaja_local
#   log(mu)     = media_global + ataque_visitante + defensa_local
#
# -----------------------------------------------------------------------------
# IDENTIFICACIÓN DEL MODELO
# -----------------------------------------------------------------------------
# La parametrización tiene DOS direcciones planas, no una. Escribiendo
# m = media_global, α = ataques, β = defensas, γ = ventaja local:
#
#   (i)  α_i → α_i + c ,  β_i → β_i − c        deja λ y μ invariantes
#   (ii) m   → m + d   ,  β_i → β_i − d        deja λ y μ invariantes
#
# Fijar α₁ = 0 remueve (i), pero NO remueve (ii): la transformación
# "subir m en d y bajar todas las defensas en d" respeta α₁ = 0 y da
# exactamente la misma verosimilitud. Verificado numéricamente en la Sección 4.
#
# Efecto práctico: media_global y el NIVEL de ataques y defensas son arbitrarios
# (dependen de dónde arrancó el optimizador, no de los datos). Quedan bien
# identificados: λ, μ, todas las probabilidades, la log-verosimilitud, el test
# LR sobre ρ, y las DIFERENCIAS entre equipos (o sea, los rankings).
#
# Se resuelve con una normalización POSTERIOR a la estimación (Sección 3.1),
# imponiendo DOS restricciones de escala: media(ataques) = 0 y media(defensas) = 0.
# Ninguna altera λ, μ ni la verosimilitud.
#
# Por qué dos y no una: α₁ = 0 es también una elección arbitraria, y deja a
# media_global valiendo la tasa del PRIMER EQUIPO ALFABÉTICO en lugar de la de
# un equipo promedio. Centrar los ataques la reemplaza por una referencia
# interpretable y hace que los parámetros se lean como desvíos respecto de la liga.
#
# Se optó por normalizar después en vez de reparametrizar para no tocar el
# optimizador: la estimación es la misma y solo cambia la lectura de los parámetros.
# -----------------------------------------------------------------------------
#
# El parámetro ρ (rho) corrige la dependencia entre goles en marcadores bajos:
# 0-0, 1-0, 0-1 y 1-1. Para todos los demás marcadores la corrección es 1.

# Función auxiliar: factor de corrección τ para cada marcador
tau <- function(goles_l, goles_v, lambda, mu, rho) {
  if      (goles_l == 0 & goles_v == 0) { 1 - lambda * mu * rho }
  else if (goles_l == 1 & goles_v == 0) { 1 + mu * rho          }
  else if (goles_l == 0 & goles_v == 1) { 1 + lambda * rho      }
  else if (goles_l == 1 & goles_v == 1) { 1 - rho               }
  else                                  { 1                      }
}

log_verosimilitud_dc <- function(params, datos, equipos) {

  n_equipos <- length(equipos)

  media_global  <- params[1]
  ventaja_local <- params[2]
  rho           <- params[3]

  # Ataque del primer equipo fijo en 0 (restricción de identificación (i))
  ataques  <- c(0, params[4:(n_equipos + 2)])
  names(ataques) <- equipos

  # Defensas
  defensas <- params[(n_equipos + 3):(2 * n_equipos + 2)]
  names(defensas) <- equipos

  log_lik <- 0

  for (i in seq_len(nrow(datos))) {
    local     <- datos$home_team[i]
    visitante <- datos$away_team[i]
    goles_l   <- datos$fthg[i]
    goles_v   <- datos$ftag[i]

    lambda <- exp(media_global + ataques[local]     + defensas[visitante] + ventaja_local)
    mu     <- exp(media_global + ataques[visitante] + defensas[local])

    # Corrección τ de Dixon-Coles para marcadores bajos
    correccion <- tau(goles_l, goles_v, lambda, mu, rho)

    # Dominio inválido de τ: se penaliza con una función CRECIENTE en la
    # violación, no con una constante. Con una constante el gradiente numérico
    # ahí es cero y BFGS puede quedar clavado reportando convergencia.
    # Con estos datos nunca se dispara (ρ ≈ 0), pero blinda el script.
    if (correccion <= 0) return(1e6 + 1e6 * abs(correccion))

    log_lik <- log_lik +
      log(correccion) +
      dpois(goles_l, lambda, log = TRUE) +
      dpois(goles_v, mu,     log = TRUE)
  }

  return(-log_lik)  # optim() minimiza, por eso devolvemos el negativo
}


# -----------------------------------------------------------------------------
# 3. OPTIMIZACIÓN
# -----------------------------------------------------------------------------

params_init_dc <- c(
  log(mean(datos_train$fthg)),  # media global
  0.1,                          # ventaja local
  -0.1,                         # rho (Dixon & Coles encontraron valor negativo)
  rep(0, n_equipos - 1),        # ataques (n-1 porque el 1ro está fijo en 0)
  rep(0, n_equipos)             # defensas
)

cat("\nOptimizando (puede tardar 2-3 minutos)...\n")

resultado_dc <- optim(
  par     = params_init_dc,
  fn      = log_verosimilitud_dc,
  datos   = datos_train,
  equipos = equipos,
  method  = "BFGS",
  control = list(maxit = 5000, trace = 1)
)

cat("\n¿Convergió?", if (resultado_dc$convergence == 0) "SÍ ✓" else "NO - revisar", "\n")
cat("Log-verosimilitud:", round(-resultado_dc$value, 2), "\n")

# Extraer parámetros (crudos, sin normalizar)
params_dc     <- resultado_dc$par
media_global  <- params_dc[1]
ventaja_local <- params_dc[2]
rho           <- params_dc[3]
ataques       <- c(0, params_dc[4:(n_equipos + 2)])
defensas      <- params_dc[(n_equipos + 3):(2 * n_equipos + 2)]
names(ataques)  <- equipos
names(defensas) <- equipos


# -----------------------------------------------------------------------------
# 3.1 NORMALIZACIÓN DE IDENTIFICACIÓN
# -----------------------------------------------------------------------------
# La parametrización tiene dos grados de libertad no identificados (ver Sección 2).
# Se los fija con dos restricciones de escala, ambas convencionales:
#
#     media(ataques) = 0     y     media(defensas) = 0
#
# Así media_global pasa a ser el nivel de referencia de un equipo PROMEDIO
# (no el del primer equipo alfabético, que es lo que impone α₁ = 0 sin normalizar),
# y ataques y defensas se leen como desvíos respecto de la media de la liga:
#   ataque  > 0 → marca más que el promedio
#   defensa < 0 → recibe menos que el promedio
#
# λ, μ, la verosimilitud y ρ quedan EXACTAMENTE iguales.

c_ataque  <- mean(ataques)
c_defensa <- mean(defensas)

ataques      <- ataques  - c_ataque
defensas     <- defensas - c_defensa
media_global <- media_global + c_ataque + c_defensa

cat("\n--- Normalización de identificación ---\n")
cat("Desplazamiento de ataques:  ", round(c_ataque,  4), "\n")
cat("Desplazamiento de defensas: ", round(c_defensa, 4), "\n")
cat("Medias tras normalizar → ataques:", round(mean(ataques), 10),
    "| defensas:", round(mean(defensas), 10), "(ambas deben ser 0)\n")

cat("\n--- Parámetros globales (normalizados) ---\n")
cat("Media global (log):", round(media_global, 4), "\n")
cat("  → equipo promedio de visitante:", round(exp(media_global), 3), "goles\n")
cat("  → equipo promedio de local:    ",
    round(exp(media_global + ventaja_local), 3), "goles\n")
cat("Ventaja local (log):", round(ventaja_local, 4),
    "→ factor multiplicativo:", round(exp(ventaja_local), 3), "\n")
cat("Rho (corrección τ):", round(rho, 4), "\n")

# --- Contraste: medias AJUSTADAS vs observadas -------------------------------
# Se comparan los promedios de λ y μ sobre los 1140 partidos contra los goles
# observados. Con link logarítmico e intercepto, las ecuaciones de score fuerzan
# la coincidencia casi exacta. Una discrepancia apreciable indicaría un error
# de indexación en la extracción de parámetros, no un problema de ajuste.
#
# NOTA IMPORTANTE PARA LA REDACCIÓN
# No confundir esta cantidad con exp(media_global + ventaja_local), que es la
# tasa del par PROMEDIO contra PROMEDIO. Esa última queda sistemáticamente POR
# DEBAJO del promedio observado, por desigualdad de Jensen: la exponencial es
# convexa y los equipos son heterogéneos, de modo que
#
#     media(exp(α_i + β_j))  ≥  exp(media(α_i + β_j))  =  exp(0) = 1
#
# La brecha (del orden del 6% con estos datos) NO es desajuste del modelo: es
# la heterogeneidad de calidad entre equipos. Son dos cantidades distintas y
# ambas correctas. Reportarlas como si fueran comparables sería un error.

lambdas_aj <- exp(media_global + ataques[datos_train$home_team] +
                  defensas[datos_train$away_team] + ventaja_local)
mus_aj     <- exp(media_global + ataques[datos_train$away_team] +
                  defensas[datos_train$home_team])

cat("\n--- Contraste: medias ajustadas vs observadas ---\n")
cat(sprintf("  Goles de local     - ajustado: %.4f | observado: %.4f | dif: %+.4f\n",
            mean(lambdas_aj), mean(datos_train$fthg),
            mean(lambdas_aj) - mean(datos_train$fthg)))
cat(sprintf("  Goles de visitante - ajustado: %.4f | observado: %.4f | dif: %+.4f\n",
            mean(mus_aj), mean(datos_train$ftag),
            mean(mus_aj) - mean(datos_train$ftag)))

if (max(abs(c(mean(lambdas_aj) - mean(datos_train$fthg),
              mean(mus_aj)     - mean(datos_train$ftag)))) > 0.02) {
  warning("Las medias ajustadas no reproducen las observadas. ",
          "Revisar los índices de extracción de ataques y defensas.")
}

# Referencia del par promedio contra promedio (cantidad DISTINTA, ver nota)
cat(sprintf("\n  Par promedio-promedio → local: %.3f | visitante: %.3f\n",
            exp(media_global + ventaja_local), exp(media_global)))
cat("  (Por debajo de las medias observadas por Jensen, no por desajuste.)\n")


# -----------------------------------------------------------------------------
# 4. CHEQUEO DE VERIFICACIÓN: la dirección plana existe
# -----------------------------------------------------------------------------
# Se evalúa la verosimilitud en el óptimo y en puntos desplazados a lo largo de
# la dirección (ii): media_global + d, todas las defensas − d, ataques intactos.
# Si la log-verosimilitud no cambia, la dirección es plana y la normalización
# de la Sección 3.1 no es cosmética: es lo que vuelve reportables los números.
#
# Cuesta 4 evaluaciones de la función (segundos), no re-optimiza nada.

idx_defensas <- (n_equipos + 3):(2 * n_equipos + 2)
nll_opt      <- log_verosimilitud_dc(params_dc, datos_train, equipos)

chequeo_plano <- map_dfr(c(-1.3, -0.5, 0.5, 2.0), function(d) {
  p        <- params_dc
  p[1]     <- p[1] + d
  p[idx_defensas] <- p[idx_defensas] - d
  tibble(d = d,
         nll_desplazada = log_verosimilitud_dc(p, datos_train, equipos),
         diferencia     = log_verosimilitud_dc(p, datos_train, equipos) - nll_opt)
})

cat("\n--- Chequeo de verificación: dirección plana de identificación ---\n")
cat("NLL en el óptimo:", format(nll_opt, digits = 12), "\n")
print(chequeo_plano)

if (max(abs(chequeo_plano$diferencia)) < 1e-6) {
  cat("Confirmado: la dirección es plana (diferencias a nivel de error numérico).\n")
  cat("-> Sin la normalización de la Sección 3.1, media_global y el nivel de\n")
  cat("   las defensas serían arbitrarios y no reproducibles entre corridas.\n")
} else {
  warning("La dirección no resultó plana. Revisar los índices de idx_defensas ",
          "antes de usar los parámetros normalizados.")
}


# -----------------------------------------------------------------------------
# 5. GUARDAR PARA EL SCRIPT 03
# -----------------------------------------------------------------------------
# Se guarda resultado_dc SIN modificar: el Script 03 usa $value (verosimilitud,
# invariante a la normalización) y $par[3] (ρ, también invariante).
# Los parámetros normalizados se guardan aparte, para reporte.

saveRDS(
  list(
    resultado_dc          = resultado_dc,
    log_verosimilitud_dc  = log_verosimilitud_dc,
    tau                   = tau,
    datos_train           = datos_train,
    datos_test            = datos_test,
    equipos               = equipos,
    n_equipos             = n_equipos,
    rho                   = rho,
    # --- parámetros normalizados: media(ataques) = 0 y media(defensas) = 0 ---
    media_global_norm     = media_global,
    ventaja_local         = ventaja_local,
    ataques_norm          = ataques,
    defensas_norm         = defensas,
    c_ataque              = c_ataque, 
    c_defensa             = c_defensa
  ),
  "resultados_dc.rds"
)

saveRDS(partidos_por_equipo, "partidos_por_equipo.rds")

cat("\nArchivos guardados:\n")
cat("  'resultados_dc.rds'         -> usar en Script 03\n")
cat("  'partidos_por_equipo.rds'   -> tabla para el marco metodológico\n")
