# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 01: Carga de múltiples temporadas
# =============================================================================
# install.packages(c("tidyverse", "janitor"))
# =============================================================================
#
# DECISIONES DE DISEÑO
#
#   (a) DOS SERIES DE CUOTAS. Se cargan las cuotas de Bet365 relevadas varios
#       días antes del encuentro (b365h, b365d, b365a) y las de cierre
#       (b365ch, b365cd, b365ca), que football-data.co.uk publica desde la
#       temporada 2019/20. Las primeras NO son cuotas de cierre: se recogen el
#       viernes a la tarde (martes para partidos entre semana). Esa distinción
#       condiciona la lectura del overround y de la eficiencia del mercado, y
#       por eso ambas series se reportan.
#
#   (b) Las cuotas de cierre NO intervienen en el filter() de NA.
#       Esto es deliberado: el filtro de faltantes se aplica sobre las
#       variables que definen la muestra de estimación (goles) y sobre la
#       serie primaria de cuotas (pre-cierre); las cuotas de cierre son un
#       insumo de robustez y su ausencia no invalida una observación.
#
#   (c) Se crea el objeto datos_mercado con las CUATRO temporadas unidas.
#       Es exclusivo del análisis de mercado del Script 07.
#
#       REGLA DE VENTANA TEMPORAL (vale para todo el trabajo):
#         - Análisis del mercado (overround, FLB)  -> cuatro temporadas.
#           Justificación: no involucra al modelo, solo cuotas y resultados
#           observados, por lo que no necesita ser fuera de muestra.
#         - Comparación modelo vs. mercado         -> solo validación 2024/25.
#           Justificación: involucra al modelo y debe ser fuera de muestra.
# =============================================================================

library(tidyverse)
library(janitor)

# -----------------------------------------------------------------------------
# 1. CARGA DE LOS ARCHIVOS
# -----------------------------------------------------------------------------

# Columnas requeridas (si falta alguna, la carga debe fallar)
columnas_base <- c(
  "date", "home_team", "away_team",
  "fthg", "ftag", "ftr",
  "b365h", "b365d", "b365a"
)

# Columnas de cierre (opcionales: si una temporada no las trae, se completan
# con NA y se deja constancia; no se interrumpe la carga)
columnas_cierre <- c("b365ch", "b365cd", "b365ca")

cargar_temporada <- function(ruta, temporada) {

  crudo <- read_csv(ruta, show_col_types = FALSE) |>
    clean_names()

  # Verificación de columnas base
  faltan_base <- setdiff(columnas_base, names(crudo))
  if (length(faltan_base) > 0) {
    stop("Temporada ", temporada, ": faltan columnas base -> ",
         paste(faltan_base, collapse = ", "))
  }

  # Verificación de columnas de cierre
  faltan_cierre <- setdiff(columnas_cierre, names(crudo))
  if (length(faltan_cierre) > 0) {
    warning("Temporada ", temporada, ": no trae las columnas de cierre -> ",
            paste(faltan_cierre, collapse = ", "),
            ". Se completan con NA.")
    crudo[faltan_cierre] <- NA_real_
  }

  crudo |>
    select(all_of(columnas_base), all_of(columnas_cierre)) |>
    mutate(date = lubridate::dmy(date)) |>  
    filter(
      !is.na(date),
      !is.na(fthg), !is.na(ftag), !is.na(ftr),
      !is.na(b365h), !is.na(b365d), !is.na(b365a)
    ) |>
    mutate(temporada = temporada)
}

# Cargar las rutas con las ubicaciones reales de los archivos
datos_2122 <- cargar_temporada("E0_2122.csv", "2021/22")
datos_2223 <- cargar_temporada("E0_2223.csv", "2022/23")
datos_2324 <- cargar_temporada("E0_2324.csv", "2023/24")
datos_2425 <- cargar_temporada("E0_2425.csv", "2024/25")


# -----------------------------------------------------------------------------
# 2. SEPARAR EN ENTRENAMIENTO Y VALIDACIÓN
# -----------------------------------------------------------------------------

# Entrenamiento: 3 temporadas históricas
datos_train <- bind_rows(datos_2122, datos_2223, datos_2324)

# Validación: temporada actual
datos_test  <- datos_2425

cat("Partidos de entrenamiento (21/22 + 22/23 + 23/24):", nrow(datos_train), "\n")
cat("Partidos de validación   (24/25):                 ", nrow(datos_test),  "\n")

# --- CHEQUEO DE VERIFICACIÓN 1 -----------------------------------------------
# La cantidad de filas debe coincidir con la esperada: 380 partidos por temporada.

n_train_esperado <- 1140   # 3 temporadas × 380 partidos
n_test_esperado  <- 380

cat("\n--- Chequeo de verificación 1: cantidad de filas ---\n")
cat(sprintf("  datos_train: %d (esperado %d) %s\n",
            nrow(datos_train), n_train_esperado,
            if (nrow(datos_train) == n_train_esperado) "OK" else "REVISAR"))
cat(sprintf("  datos_test:  %d (esperado %d) %s\n",
            nrow(datos_test), n_test_esperado,
            if (nrow(datos_test) == n_test_esperado) "OK" else "REVISAR"))

if (nrow(datos_train) != n_train_esperado || nrow(datos_test) != n_test_esperado) {
  warning("La cantidad de filas cambió respecto de lo esperado. ",
          "Revisar el filter() antes de continuar: los Scripts 02 a 06 ",
          "dejarían de reproducir sus resultados.")
}

# --- CHEQUEO DE VERIFICACIÓN 2: integridad de los datos ----------------------
datos_todas <- bind_rows(datos_2122, datos_2223, datos_2324, datos_2425)

# Cada par (local, visitante) una sola vez por temporada
dup <- datos_todas |> count(temporada, home_team, away_team) |> filter(n > 1)
stopifnot(nrow(dup) == 0)

# Coherencia entre marcador y FTR
incoherentes <- datos_todas |>
  filter(ftr != case_when(fthg > ftag ~ "H", fthg < ftag ~ "A", TRUE ~ "D"))
stopifnot(nrow(incoherentes) == 0)

# Cuotas válidas (una cuota <= 1 daría probabilidad implícita >= 1)
stopifnot(min(c(datos_todas$b365h, datos_todas$b365d, datos_todas$b365a)) > 1)

cat("\n--- Chequeo de verificación 2: integridad --- OK\n")


# -----------------------------------------------------------------------------
# 3. DISPONIBILIDAD DE LAS CUOTAS DE CIERRE
# -----------------------------------------------------------------------------
# Este recuento va al marco metodológico: hay que declarar explícitamente
# cuántos partidos no tienen cuota de cierre y en qué temporadas.

cobertura_cierre <- datos_todas |>
  summarise(
    partidos           = n(),
    con_cierre         = sum(!is.na(b365ch) & !is.na(b365cd) & !is.na(b365ca)),
    sin_cierre         = sum( is.na(b365ch) |  is.na(b365cd) |  is.na(b365ca)),
    .by = temporada
  ) |>
  mutate(cobertura_pct = round(con_cierre / partidos * 100, 2)) |>
  arrange(temporada)

cat("\n--- Cobertura de las cuotas de cierre (B365CH/CD/CA) ---\n")
print(cobertura_cierre)

cat(sprintf(
  "\nTotal: %d partidos, %d con cuota de cierre (%.2f%%), %d sin cuota de cierre.\n",
  nrow(datos_todas),
  sum(cobertura_cierre$con_cierre),
  sum(cobertura_cierre$con_cierre) / nrow(datos_todas) * 100,
  sum(cobertura_cierre$sin_cierre)
))
cat("(Este número va al marco metodológico.)\n")


# -----------------------------------------------------------------------------
# 4. VERIFICAR QUE LOS EQUIPOS DE 24/25 ESTÉN EN EL ENTRENAMIENTO
# -----------------------------------------------------------------------------
# Algunos equipos pueden haber ascendido entre temporadas.
# El modelo solo puede predecir equipos que "conoce" del entrenamiento.

equipos_train <- unique(c(datos_train$home_team, datos_train$away_team))
equipos_test  <- unique(c(datos_test$home_team,  datos_test$away_team))

equipos_nuevos <- setdiff(equipos_test, equipos_train)

if (length(equipos_nuevos) == 0) {
  cat("\nTodos los equipos de 24/25 están en el entrenamiento ✓\n")
} else {
  cat("\nEquipos en 24/25 que NO estaban en temporadas anteriores:\n")
  print(equipos_nuevos)
  cat("-> Estos equipos ascendieron. Los vamos a manejar en el Paso 2.\n")
}


# -----------------------------------------------------------------------------
# 5. ANÁLISIS EXPLORATORIO COMPARATIVO ENTRE TEMPORADAS
# -----------------------------------------------------------------------------

cat("\n--- Distribución de resultados por temporada ---\n")
datos_todas |>
  count(temporada, ftr) |>
  mutate(
    resultado  = case_when(ftr == "H" ~ "Local",
                           ftr == "D" ~ "Empate",
                           ftr == "A" ~ "Visitante"),
    proporcion = round(n / sum(n), 3),
    .by = temporada
  ) |>
  select(temporada, resultado, n, proporcion) |>
  arrange(temporada, resultado) |>
  print()

cat("\n--- Media de goles por temporada ---\n")
datos_todas |>
  summarise(
    media_local       = round(mean(fthg), 3),
    media_visitante   = round(mean(ftag), 3),
    total_por_partido = round(mean(fthg + ftag), 3),
    .by = temporada
  ) |>
  print()


# -----------------------------------------------------------------------------
# 6. GRÁFICO: DISTRIBUCIÓN DE GOLES POR TEMPORADA
# -----------------------------------------------------------------------------

datos_todas |>
  pivot_longer(cols = c(fthg, ftag),
               names_to  = "equipo",
               values_to = "goles") |>
  mutate(equipo = if_else(equipo == "fthg", "Local", "Visitante")) |>
  ggplot(aes(x = goles, fill = equipo)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7, color = "white") +
  facet_wrap(~temporada, nrow = 2) +
  scale_fill_manual(values = c("Local" = "#1f77b4", "Visitante" = "#ff7f0e")) +
  labs(
    title    = "Distribución de goles por temporada - Premier League",
    subtitle = "Entrenamiento: 21/22 a 23/24  |  Validación: 24/25",
    x        = "Goles por partido",
    y        = "Frecuencia",
    fill     = "Equipo"
  ) +
  theme_minimal()


# -----------------------------------------------------------------------------
# 7. OBJETO PARA EL ANÁLISIS DE MERCADO
# -----------------------------------------------------------------------------
# datos_mercado reúne las cuatro temporadas. Es el insumo EXCLUSIVO del
# análisis de mercado del Script 07 (overround y favourite-longshot bias).
#
# No reemplaza a datos_train ni a datos_test: esos objetos siguen definiendo
# la partición del modelo y no se tocan.

datos_mercado <- datos_todas

cat("\n--- Objeto datos_mercado (análisis de mercado, Script 07) ---\n")
cat("Partidos:", nrow(datos_mercado), "sobre", n_distinct(datos_mercado$temporada), "temporadas\n")
cat("Regla: análisis del mercado = 4 temporadas | comparación modelo-mercado = solo 24/25\n")


# -----------------------------------------------------------------------------
# 8. GUARDAR
# -----------------------------------------------------------------------------

saveRDS(datos_train,      "datos_train.rds")
saveRDS(datos_test,       "datos_test.rds")
saveRDS(datos_mercado,    "datos_mercado.rds")
saveRDS(cobertura_cierre, "cobertura_cuotas_cierre.rds")
saveRDS(equipos_nuevos, "equipos_nuevos.rds")

cat("\nArchivos guardados:\n")
cat("  'datos_train.rds'              -> para entrenar el modelo (Paso 2)\n")
cat("  'datos_test.rds'               -> para validar el modelo  (Pasos 3 y 4)\n")
cat("  'datos_mercado.rds'            -> análisis de mercado, 4 temporadas (Script 07)\n")
cat("  'cobertura_cuotas_cierre.rds'  -> tabla de faltantes para el marco metodológico\n")
cat("  'equipos_nuevos.rds'           -> equipos que no estan en train y si en test\n")