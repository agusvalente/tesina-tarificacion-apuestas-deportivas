# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 04: Modelo de Maher (modelo definitivo)
# =============================================================================
# Depende de: Script 01 (datos_train.rds, datos_test.rds, equipos_nuevos.rds)
#             Script 03 (resultados_lrt.rds)
#
# El test de razón de verosimilitudes del Script 03 mostró que ρ no es
# estadísticamente significativo (LR = 0,053; p = 0,817; IC95 de ρ =
# (−0,066 ; 0,084)), por lo que se adopta el modelo de Maher (1982) como
# especificación final.
# =============================================================================
#
# DECISIONES DE DISEÑO
#
#   (a) SE USAN LOS PARÁMETROS NORMALIZADOS DEL SCRIPT 03 (Sección 2).
#       Los parámetros crudos de resultado_maher$par dan media_global = 0,652,
#       que NO es la media de la liga sino la tasa del primer equipo
#       alfabético (Arsenal), porque el optimizador impone α₁ = 0 como recurso
#       de identificación.
#
#       Importa para los gráficos, que dibujan líneas de referencia en cero.
#       Con los parámetros crudos ese cero marcaría a Arsenal (el tercer mejor
#       ataque de la liga) y no al promedio, y el gráfico de cuadrantes
#       sugeriría que 22 de 26 equipos están "por debajo de la referencia". Con
#       los parámetros normalizados (media(ataques) = 0 y media(defensas) = 0)
#       el cero significa lo que los gráficos dicen que significa, y el mismo
#       criterio rige para ataque y para defensa.
#
#       Las probabilidades NO dependen de esta elección: λ y μ son invariantes
#       a la normalización.
#
#   (b) LA GRILLA DE MARCADORES SE TRUNCA EN 25 GOLES (Sección 5).
#       Con una grilla corta las probabilidades no suman 1 (con max_goles = 9,
#       mínimo observado 0,998). Lo relevante no es la magnitud sino la
#       DIRECCIÓN: la masa que se pierde sale casi íntegramente del lado del
#       equipo con λ más alto. Para λ = 3,16 y μ = 0,90 el déficit de 1,6e-03
#       recae entero sobre prob_H, dejando prob_D y prob_A intactas.
#
#       O sea, una truncación corta SUBESTIMA SISTEMÁTICAMENTE la probabilidad
#       de que gane el favorito, y tanto más cuanto más marcado sea el
#       favoritismo. En una tesina sobre favourite-longshot bias eso
#       introduciría por vía computacional un sesgo con la misma forma que el
#       fenómeno estudiado.
#
#       Déficit de masa según la grilla (verificado numéricamente):
#           λ, μ         K = 9      K = 15     K = 25
#           1,55 / 1,25  6,3e-06    1,3e-11    0
#           2,40 / 1,00  2,0e-04    6,1e-09    0
#           3,16 / 0,90  1,6e-03    2,5e-07    1,6e-15
#           3,50 / 1,20  3,3e-03    9,2e-07    1,2e-14
#
#       Con K = 25 el error queda bajo la precisión de máquina, de modo que la
#       verificación de la Sección 5 es un chequeo que efectivamente puede
#       fallar.
#
#   (c) ANÁLISIS DE SENSIBILIDAD AL PESO DE IPSWICH (Sección 3.1).
#       peso_descendidos es un grado de libertad del diseño sin justificación
#       previa, e Ipswich interviene en 38 de los 380 partidos de validación
#       (el 10% de la muestra sobre la que se evalúa todo el modelo). Entre
#       peso 1 y peso 3 el ataque imputado cambia 0,104 en escala logarítmica,
#       un 11% en la tasa de goles.
#
#       El criterio no se cambia después de ver su efecto (eso sería post-hoc):
#       se documenta ese efecto. La tabla va al marco metodológico.
#
#   (d) Los equipos ascendidos no están escritos a mano en el código: se leen
#       de equipos_nuevos.rds del Script 01, con chequeo cruzado contra los
#       datos.
#
#   (e) SE DOCUMENTA EL CAMBIO DE RÉGIMEN DE LA VENTAJA LOCAL (Sección 5.1).
#       El modelo estima exp(γ) = 1,236 a partir del entrenamiento, pero en la
#       temporada de validación el cociente observado fue 1,063: la ventaja
#       local más débil de las cuatro temporadas y por margen amplio.
#       Consecuencia: el modelo sobreestima P(H) de forma sistemática en los
#       380 partidos de validación. Hay que declararlo antes de interpretar la
#       calibración del Script 08, o ese resultado se leerá como un defecto de
#       la especificación Poisson cuando en realidad es un cambio de régimen
#       entre muestras. La tabla se guarda para el marco metodológico.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. CARGAR DATOS Y RESULTADOS PREVIOS
# -----------------------------------------------------------------------------

datos_train    <- readRDS("datos_train.rds")
datos_test     <- readRDS("datos_test.rds")
lrt            <- readRDS("resultados_lrt.rds")
equipos_nuevos <- readRDS("equipos_nuevos.rds")

equipos   <- sort(unique(c(datos_train$home_team, datos_train$away_team)))
n_equipos <- length(equipos)

equipos_test <- unique(c(datos_test$home_team, datos_test$away_team))

# Chequeo cruzado: el objeto guardado por el Script 01 tiene que coincidir con
# lo que se deduce de los datos cargados acá.
stopifnot(identical(sort(equipos_nuevos), sort(setdiff(equipos_test, equipos))))

cat("Modelo de Maher cargado desde el Script 03.\n")
cat("Log-verosimilitud:", round(lrt$loglik_maher, 4), "\n")
cat("¿Convergió?", if (lrt$resultado_maher$convergence == 0) "SÍ ✓" else "NO - revisar", "\n")
cat("Equipos ascendidos a la 24/25 sin historial de entrenamiento:",
    paste(equipos_nuevos, collapse = ", "), "\n")


# -----------------------------------------------------------------------------
# 2. PARÁMETROS (NORMALIZADOS EN EL SCRIPT 03)
# -----------------------------------------------------------------------------
# Convención de identificación: media(ataques) = 0 y media(defensas) = 0.
# Ambas son restricciones de escala y no alteran λ, μ ni la verosimilitud.
# El mismo criterio rige para ataque y para defensa, de modo que el cero
# significa "promedio de la liga" en las dos escalas.
#
# Ver la Sección 2 del Script 02 para la derivación de por qué hacen falta DOS
# restricciones y no una.

media_global  <- lrt$media_global_norm
ventaja_local <- lrt$ventaja_local
ataques       <- lrt$ataques_norm
defensas      <- lrt$defensas_norm

stopifnot(abs(mean(ataques))  < 1e-8,
          abs(mean(defensas)) < 1e-8,
          identical(names(ataques),  equipos),
          identical(names(defensas), equipos))

cat("\n--- Parámetros globales (normalizados) ---\n")
cat("Media global (log):", round(media_global, 4), "\n")
cat("  → equipo promedio de local:    ",
    round(exp(media_global + ventaja_local), 3), "goles\n")
cat("  → equipo promedio de visitante:",
    round(exp(media_global), 3), "goles\n")
cat("Ventaja local (log):", round(ventaja_local, 4),
    "→ factor multiplicativo:", round(exp(ventaja_local), 3), "\n")
cat("Referencia: ataque 0 y defensa 0 = promedio de la liga (mismo criterio).\n")


# -----------------------------------------------------------------------------
# 3. IMPUTAR LOS EQUIPOS ASCENDIDOS
# -----------------------------------------------------------------------------
# Los equipos ascendidos no jugaron en las temporadas de entrenamiento, así que
# el modelo no tiene parámetros estimados para ellos.
#
# Criterio: promedio ponderado de la liga, dándole más peso a los equipos que
# participaron de la Premier entre 21/22 y 23/24 pero no están en la 24/25.
# La lógica es que esos equipos son más similares a un recién ascendido: son
# los que no lograron mantenerse en la categoría.
#
# LIMITACIÓN DEL CRITERIO (para la discusión, no para cambiar el código):
# Leicester y Southampton ascendieron a la 24/25 junto con Ipswich y reciben
# peso 1, aunque son los dos equipos más comparables a Ipswich que existen en
# la muestra. El criterio elegido es razonable pero no es el único razonable;
# la Sección 3.1 documenta cuánto depende el resultado de esa elección.
#
# PARÁMETRO DE DISEÑO:
peso_descendidos <- 3   # 1 = promedio simple | 2 = doble peso | 3 = triple, etc.

equipos_ausentes <- setdiff(equipos, equipos_test)

cat("\n--- Equipos ausentes en 24/25 (peso:", peso_descendidos, ") ---\n")
print(equipos_ausentes)

pesos <- ifelse(equipos %in% equipos_ausentes, peso_descendidos, 1)
names(pesos) <- equipos

ataque_ponderado  <- weighted.mean(ataques,  pesos)
defensa_ponderado <- weighted.mean(defensas, pesos)

# Asignación a todos los ascendidos (no hardcodeada a un nombre)
ataques[equipos_nuevos]  <- ataque_ponderado
defensas[equipos_nuevos] <- defensa_ponderado

cat("\n--- Parámetros imputados a los ascendidos ---\n")
cat("Ataque ponderado: ", round(ataque_ponderado,  4), "\n")
cat("Defensa ponderada:", round(defensa_ponderado, 4), "\n")
cat("Referencia, promedio simple de la liga → ataque:",
    round(mean(ataques[equipos]), 4),
    "| defensa:", round(mean(defensas[equipos]), 4),
    "(ambos 0 por normalización)\n")

n_partidos_nuevos <- sum(datos_test$home_team %in% equipos_nuevos |
                         datos_test$away_team %in% equipos_nuevos)
cat("Partidos de validación que involucran a un ascendido:",
    n_partidos_nuevos, "de", nrow(datos_test),
    sprintf("(%.1f%%)\n", n_partidos_nuevos / nrow(datos_test) * 100))


# -----------------------------------------------------------------------------
# 3.1 SENSIBILIDAD AL PESO DE LOS EQUIPOS AUSENTES
# -----------------------------------------------------------------------------
# peso_descendidos no tiene justificación externa: es una elección del diseño.
# Como los ascendidos intervienen en ~10% de la muestra de validación, hay que
# mostrar explícitamente cuánto mueve esa elección. Esta tabla va al marco
# metodológico.
#
# Las tasas se calculan contra un rival promedio (ataque 0, defensa 0), que con
# la normalización vigente es literalmente el equipo medio de la liga.

sensibilidad_peso <- map_dfr(c(1, 2, 3, 5, 10), function(w) {
  p <- ifelse(equipos %in% equipos_ausentes, w, 1)
  a <- weighted.mean(ataques[equipos],  p)
  d <- weighted.mean(defensas[equipos], p)
  tibble(
    peso           = w,
    ataque_imput   = a,
    defensa_imput  = d,
    goles_de_local = exp(media_global + a + 0 + ventaja_local),
    goles_de_visit = exp(media_global + a + 0)
  )
})

cat("\n--- Sensibilidad de los parámetros imputados al peso ---\n")
print(sensibilidad_peso |> mutate(across(where(is.numeric), \(x) round(x, 4))))

cat(sprintf(
  "\nRango del ataque imputado entre peso 1 y peso 10: %.4f (%.1f%% en tasa de goles).\n",
  diff(range(sensibilidad_peso$ataque_imput)),
  (exp(diff(range(sensibilidad_peso$ataque_imput))) - 1) * 100
))
cat("-> En el Script 08, recalcular Brier y log-loss EXCLUYENDO los partidos de\n")
cat("   los ascendidos. Si las conclusiones se mantienen, el peso deja de ser\n")
cat("   relevante y se declara resuelto. Si no, es una limitación a escribir.\n")


# -----------------------------------------------------------------------------
# 4. TABLA DE PARÁMETROS
# -----------------------------------------------------------------------------

equipos_todos <- c(equipos, equipos_nuevos)

tabla_params <- tibble(
  equipo   = equipos_todos,
  ataque   = round(ataques[equipos_todos],  4),
  defensa  = round(defensas[equipos_todos], 4),
  imputado = equipos_todos %in% equipos_nuevos
) |>
  arrange(desc(ataque))

cat("\n--- Parámetros por equipo (0 = promedio de la liga) ---\n")
print(tabla_params, n = Inf)


# -----------------------------------------------------------------------------
# 5. PROBABILIDADES PARA LOS PARTIDOS DE VALIDACIÓN (24/25)
# -----------------------------------------------------------------------------
# Para cada partido se suman las probabilidades de todos los marcadores de la
# grilla 0..max_goles en cada equipo.
#
# max_goles = 25. Con 9 la masa truncada llegaría a 3,3e-03 y recaería
# ENTERA sobre la probabilidad del favorito, subestimándola de forma sistemática
# — un sesgo con la misma forma que el favourite-longshot bias que la tesina
# mide. Con 25 el déficit queda bajo la precisión de máquina. Ver nota (b).

calcular_probs <- function(lambda, mu, max_goles = 25) {
  prob_matriz <- outer(
    0:max_goles, 0:max_goles,
    FUN = function(i, j) dpois(i, lambda) * dpois(j, mu)
  )
  c(
    prob_H = sum(prob_matriz[lower.tri(prob_matriz)]),  # filas = goles local
    prob_D = sum(diag(prob_matriz)),
    prob_A = sum(prob_matriz[upper.tri(prob_matriz)])   # columnas = goles visitante
  )
}

cat("\nCalculando probabilidades para los", nrow(datos_test), "partidos de 24/25...\n")

resultados <- datos_test |>
  rowwise() |>
  mutate(
    lambda = exp(media_global + ataques[home_team] + defensas[away_team] + ventaja_local),
    mu     = exp(media_global + ataques[away_team] + defensas[home_team]),
    probs  = list(calcular_probs(lambda, mu)),
    prob_H = probs[["prob_H"]],
    prob_D = probs[["prob_D"]],
    prob_A = probs[["prob_A"]]
  ) |>
  select(-probs) |>
  ungroup()

cat("\nPrimeros 6 partidos con probabilidades del modelo:\n")
resultados |>
  select(home_team, away_team, ftr, lambda, mu, prob_H, prob_D, prob_A) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  head(6) |>
  print()

# --- CHEQUEO DE VERIFICACIÓN: las probabilidades tienen que sumar 1 -----------
# Con max_goles = 25 este chequeo puede FALLAR de verdad. Con max_goles = 9
# pasaría con un desvío de hasta 2e-03 concentrado en la probabilidad del favorito.
suma_probs   <- resultados$prob_H + resultados$prob_D + resultados$prob_A
desvio_max   <- max(abs(suma_probs - 1))

cat("\n--- Chequeo de verificación: suma de probabilidades ---\n")
cat("Desvío máximo respecto de 1:", format(desvio_max, digits = 3), "\n")

if (desvio_max > 1e-10) {
  warning("Las probabilidades del modelo no suman 1 (desvío ", desvio_max, "). ",
          "Aumentar max_goles antes de usarlas en los Scripts 08 y 09: las ",
          "proper scoring rules exigen vectores de probabilidad completos.")
} else {
  cat("Truncación de la grilla por debajo de la precisión numérica ✓\n")
}

# --- Chequeo de coherencia: λ y μ agregados vs goles observados ---------------
# Es fuera de muestra, así que no tienen por qué coincidir exactamente; sirve
# para detectar errores de indexación, no para evaluar ajuste.
cat(sprintf(
  "\nGoles esperados vs observados en validación → local: %.3f / %.3f | visitante: %.3f / %.3f\n",
  mean(resultados$lambda), mean(datos_test$fthg),
  mean(resultados$mu),     mean(datos_test$ftag)))


# -----------------------------------------------------------------------------
# 5.1 VENTAJA LOCAL: ENTRENAMIENTO VS VALIDACIÓN
# -----------------------------------------------------------------------------
# El chequeo anterior no falla, pero sus dos desvíos van en direcciones opuestas
# y con magnitud casi idéntica: el modelo sobreestima los goles de local y
# subestima los de visitante en proporciones similares. Eso no es ruido, es la
# ventaja local.
#
# El modelo la estima como una CONSTANTE a partir del entrenamiento. Si la
# temporada de validación operó en otro régimen, el modelo va a errar de forma
# sistemática y direccional, sobreestimando P(H) en todos los partidos.
#
# POR QUÉ IMPORTA (marco metodológico y discusión de resultados):
#
#   1. En el Script 08 esto va a aparecer como calibración deficiente del
#      modelo en la categoría local. Leerlo como "el Poisson no calibra bien"
#      sería incorrecto: la causa es un cambio de régimen entre las muestras,
#      no un defecto de la especificación.
#
#   2. Da la explicación mecánica de por qué el mercado le gana al modelo en
#      RESOLUCIÓN: el mercado incorpora la ventaja local vigente partido a
#      partido, mientras el modelo la tiene congelada en el promedio de tres
#      temporadas anteriores.
#
#   3. Anticipa la objeción obvia sobre la ponderación temporal φ(t) de Dixon &
#      Coles que este trabajo no implementa. Mirando la serie de cocientes por
#      temporada NO hay tendencia: 24/25 es un quiebre, no la continuación de
#      una caída. Un modelo con decaimiento exponencial habría pesado más la
#      23/24 y tampoco habría anticipado el valor de 24/25. La omisión de φ(t)
#      sigue siendo una limitación, pero NO es la causa de este desvío.
#
# ADVERTENCIA: no re-estimar la ventaja local con datos de 24/25. Sería
# look-ahead bias y rompería el diseño fuera de muestra. Esto se documenta y se
# usa para interpretar, no para corregir el modelo.

ventaja_por_temporada <- bind_rows(
  datos_train |> mutate(muestra = "Entrenamiento"),
  datos_test  |> mutate(muestra = "Validación")
) |>
  summarise(
    partidos        = n(),
    goles_local     = mean(fthg),
    goles_visitante = mean(ftag),
    cociente        = mean(fthg) / mean(ftag),
    pct_victoria_h  = mean(ftr == "H"),
    .by = c(muestra, temporada)
  ) |>
  arrange(temporada)

cat("\n--- Ventaja local por temporada ---\n")
print(ventaja_por_temporada |> mutate(across(where(is.double), \(x) round(x, 3))))

vent_modelo <- exp(ventaja_local)
vent_train  <- mean(datos_train$fthg) / mean(datos_train$ftag)
vent_test   <- mean(datos_test$fthg)  / mean(datos_test$ftag)
brecha_pct  <- (vent_test / vent_modelo - 1) * 100

ventaja_local_regimen <- tibble(
  fuente    = c("Estimada por el modelo (21/22-23/24)",
                "Observada en entrenamiento",
                "Observada en validación (24/25)"),
  cociente  = c(vent_modelo, vent_train, vent_test)
) |>
  mutate(brecha_vs_modelo_pct = (cociente / vent_modelo - 1) * 100)

cat("\n--- Ventaja local: entrenamiento vs validación ---\n")
print(ventaja_local_regimen |> mutate(across(where(is.double), \(x) round(x, 4))))

cat(sprintf("\nBrecha del modelo respecto de la validación: %+.1f%%\n", brecha_pct))

if (abs(brecha_pct) > 5) {
  cat("\nATENCIÓN: la ventaja local cambió de régimen entre entrenamiento y\n")
  cat("validación. El modelo sobreestimará P(H) de forma sistemática en los 380\n")
  cat("partidos. Es una limitación del modelo estático, no un defecto de la\n")
  cat("especificación Poisson. DECLARARLO antes de interpretar la calibración\n")
  cat("del Script 08 y la descomposición de Murphy del Script 09.\n")
} else {
  cat("La ventaja local se mantiene estable entre muestras ✓\n")
}


# -----------------------------------------------------------------------------
# 6. GRÁFICOS
# -----------------------------------------------------------------------------
# Con los parámetros normalizados, la línea en cero marca el promedio de la liga
# tanto en ataque como en defensa.

etiqueta_imput <- paste0("Naranja = imputado (",
                         paste(equipos_nuevos, collapse = ", "), ")")
subtitulo <- paste0("Modelo Maher | Entrenado con 21/22 - 23/24 | 0 = promedio de la liga | ",
                    etiqueta_imput)

# --- Ranking de ataque ---
tabla_params |>
  ggplot(aes(x = reorder(equipo, ataque), y = ataque, fill = imputado)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("FALSE" = "#e74c3c", "TRUE" = "#f39c12"),
                    guide  = "none") +
  coord_flip() +
  labs(
    title    = "Parámetro de ataque por equipo",
    subtitle = subtitulo,
    x        = NULL,
    y        = "Parámetro de ataque (0 = promedio de la liga; mayor = más goleador)"
  ) +
  theme_minimal()

# --- Ranking de defensa ---
tabla_params |>
  ggplot(aes(x = reorder(equipo, -defensa), y = defensa, fill = imputado)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("FALSE" = "#2980b9", "TRUE" = "#f39c12"),
                    guide  = "none") +
  coord_flip() +
  labs(
    title    = "Parámetro de defensa por equipo",
    subtitle = subtitulo,
    x        = NULL,
    y        = "Parámetro de defensa (0 = promedio de la liga; más negativo = mejor)"
  ) +
  theme_minimal()

# --- Mapa ataque vs defensa ---
tabla_params |>
  ggplot(aes(x = defensa, y = ataque, label = equipo, color = imputado)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#f39c12"),
                     guide  = "none") +
  labs(
    title    = "Ataque vs defensa por equipo",
    subtitle = subtitulo,
    x        = "Parámetro de defensa (0 = promedio; más negativo = mejor)",
    y        = "Parámetro de ataque (0 = promedio; más positivo = mejor)"
  ) +
  theme_minimal()


# -----------------------------------------------------------------------------
# 7. GUARDAR PARA LOS SCRIPTS SIGUIENTES
# -----------------------------------------------------------------------------
# Se conservan los nombres de los elementos que ya leen los Scripts 05 a 09.

saveRDS(
  list(
    ataques           = ataques,
    defensas          = defensas,
    media_global      = media_global,
    ventaja_local     = ventaja_local,
    equipos           = equipos,
    tabla_params      = tabla_params,
    resultados        = resultados,
    # --- agregados ---
    equipos_nuevos        = equipos_nuevos,
    equipos_ausentes      = equipos_ausentes,
    peso_descendidos      = peso_descendidos,
    sensibilidad_peso     = sensibilidad_peso,
    max_goles             = 25,
    # --- régimen de ventaja local (interpretación de los Scripts 08 y 09) ---
    ventaja_local_regimen = ventaja_local_regimen,
    ventaja_por_temporada = ventaja_por_temporada,
    vent_modelo           = vent_modelo,
    vent_test             = vent_test,
    brecha_ventaja_pct    = brecha_pct
  ),
  "resultados_maher.rds"
)

saveRDS(ventaja_por_temporada, "ventaja_local_por_temporada.rds")

cat("\nArchivos guardados:\n")
cat("  'resultados_maher.rds'              -> usar en Scripts 05 y 06\n")
cat("  'ventaja_local_por_temporada.rds'   -> tabla para el marco metodológico\n")

cat("\nPendientes anotados para los scripts siguientes:\n")
cat("  Script 08: recalcular Brier y log-loss EXCLUYENDO los partidos de los\n")
cat("             ascendidos (sensibilidad al peso, Sección 3.1).\n")
cat("  Script 08: declarar el cambio de régimen de la ventaja local (Sección 5.1)\n")
cat("             ANTES de interpretar la calibración del modelo en la categoría H.\n")
cat("  Script 09: usar la brecha de ventaja local para explicar la diferencia de\n")
cat("             resolución entre modelo y mercado.\n")
