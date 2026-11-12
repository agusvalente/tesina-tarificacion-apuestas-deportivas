# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 09: Síntesis final
# =============================================================================
# install.packages(c("tidyverse", "ggrepel"))
# =============================================================================
# Depende de: Script 07 (resultados_mercado.rds)
#             Script 08 (resultados_calibracion.rds)
#
# Objetivos:
#   1. Construir la tabla maestra de resultados y exportarla como .csv
#   2. Regenerar y guardar todos los gráficos como .png publicables (300 dpi)
#   3. Producir el gráfico de síntesis: modelo vs. mercado partido a partido
#   4. Identificar y exportar los 10 partidos con mayor y menor discrepancia
#      entre las probabilidades del modelo y del mercado
#   5. Tabla y gráfico de parámetros del modelo final (ataque y defensa)
#   6. Top 10 partidos donde el modelo acertó más y falló más
#   7. Estadísticas descriptivas de la temporada 24/25
#   8. Exportar los CONTRASTES PRINCIPALES de los Scripts 07 y 08
#      (Sección 2.1), sin los cuales la síntesis entregaría magnitudes sin
#      inferencia
#
# -----------------------------------------------------------------------------
# DECISIONES DE DISEÑO
#
#   (a) LÍNEA DE REFERENCIA DEL GRÁFICO DE CAIN (07-F). No está en cero sino
#       en el retorno esperado bajo calibración perfecta y margen uniforme,
#       1/(1+OR) − 1 (mercado$referencia_justa), igual que en el Script 07.
#       Con overround NINGUNA apuesta puede rendir cero, de modo que una línea
#       en cero haría parecer que las cuatro bandas pierden dinero cuando lo
#       que corresponde juzgar es su desvío respecto de ese retorno esperado.
#
#   (b) ETIQUETADO DEL GRÁFICO DE PARÁMETROS (09-B). Los parámetros se estiman
#       sobre la ventana de entrenamiento 2021/22–2023/24 y la tabla incluye
#       26 equipos, seis de los cuales NO disputaron la 24/25 (descendidos en
#       temporadas previas). El título declara esa ventana, y las etiquetas
#       distinguen a los equipos no presentes en 24/25 y a Ipswich, cuyos
#       parámetros no están estimados sino IMPUTADOS vía peso_descendidos.
#
#   (c) CONTRASTES PRINCIPALES (Sección 2.1). Junto con las desagregaciones
#       se exportan los contrastes que los Scripts 07 y 08 declaran
#       instrumento de referencia: el test agregado de FLB con errores
#       estándar agrupados, el test de tendencia económica, la robustez
#       temporal del pareado de overround, el contraste primario de
#       calibración y el piso de ruido de REL.
#
#       El caso crítico es el FLB: tabla_flb_regresiones.csv incluye una fila
#       que declara "b > 1: consistente con FLB" —Visitante en la serie de
#       cierre, p = 0,045— y el Script 07 advierte expresamente que NO hay que
#       encabezar con ese rechazo aislado. Exportar las seis regresiones sin
#       el test agregado al lado llevaría al lector exactamente a la
#       conclusión que el Script 07 quiere evitar.
#
#   (d) ADVERTENCIAS HEREDADAS DEL SCRIPT 08. El gráfico 08-B (log-loss por
#       decil) advierte que al condicionar sobre el resultado deja de ser una
#       regla de puntuación propia, y la fila REL de la tabla maestra advierte
#       que las calibraciones son indistinguibles del piso de ruido de
#       binning.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. CARGAR DATOS
# -----------------------------------------------------------------------------

mercado      <- readRDS("resultados_mercado.rds")
calibracion  <- readRDS("resultados_calibracion.rds")
modelo_final <- readRDS("resultados_modelo_final.rds")

datos        <- mercado$datos          # 380 partidos con overround y pi_*
datos_largo  <- mercado$datos_largo    # 1.140 filas resultado × partido
datos_bins   <- mercado$datos_bins     # 10 deciles del Script 07

metricas_globales      <- calibracion$metricas_globales
descomposicion         <- calibracion$descomposicion
metricas_por_resultado <- calibracion$metricas_por_resultado
metricas_por_bin       <- calibracion$metricas_por_bin

cat("Datos cargados correctamente.\n")
cat("Partidos:", nrow(datos), "| Observaciones resultado-partido:", nrow(datos_largo), "\n")


# -----------------------------------------------------------------------------
# 2. TABLA MAESTRA DE RESULTADOS
# -----------------------------------------------------------------------------
# Consolida los números clave de la comparación modelo vs. mercado.
# Una fila por métrica, una columna por competidor.
#
# ATENCIÓN A LA VENTANA TEMPORAL: esta tabla es de la comparación modelo vs.
# mercado y por lo tanto corresponde a la validación 2024/25 únicamente.
# El overround que figura acá es el de 24/25. El overround del ANÁLISIS DE
# MERCADO (cuatro temporadas) va en la tabla de la Sección 2b.

leer <- function(tabla, columna, fuente) {
  tabla[[columna]][tabla$fuente == fuente]
}

tabla_maestra <- tibble(
  Metrica = c(
    "Overround promedio 24/25 (%)",
    "Brier Score global",
    "Brier Skill Score",
    "Log-loss global",
    "Calibración BS — REL (Murphy, 1973)",
    "Resolución BS — RES (Murphy, 1973)",
    "Incertidumbre BS — UNC (Murphy, 1973)"
  ),
  `Mercado (pre-cierre)` = c(
    round(mean(datos$overround) * 100, 4),
    round(leer(metricas_globales, "brier_score", "Mercado (pre-cierre)"), 6),
    round(leer(metricas_globales, "brier_skill", "Mercado (pre-cierre)"), 4),
    round(leer(metricas_globales, "log_loss",    "Mercado (pre-cierre)"), 6),
    round(leer(descomposicion, "calibracion",   "Mercado (pre-cierre)"),  6),
    round(leer(descomposicion, "resolucion",    "Mercado (pre-cierre)"),  6),
    round(leer(descomposicion, "incertidumbre", "Mercado (pre-cierre)"),  6)
  ),
  `Mercado (cierre)` = c(
    round(mean(datos$overround_c, na.rm = TRUE) * 100, 4),
    round(leer(metricas_globales, "brier_score", "Mercado (cierre)"), 6),
    round(leer(metricas_globales, "brier_skill", "Mercado (cierre)"), 4),
    round(leer(metricas_globales, "log_loss",    "Mercado (cierre)"), 6),
    round(leer(descomposicion, "calibracion",   "Mercado (cierre)"),  6),
    round(leer(descomposicion, "resolucion",    "Mercado (cierre)"),  6),
    round(leer(descomposicion, "incertidumbre", "Mercado (cierre)"),  6)
  ),
  Modelo = c(
    NA,   # el overround no aplica al modelo
    round(leer(metricas_globales, "brier_score", "Modelo"), 6),
    round(leer(metricas_globales, "brier_skill", "Modelo"), 4),
    round(leer(metricas_globales, "log_loss",    "Modelo"), 6),
    round(leer(descomposicion, "calibracion",   "Modelo"),  6),
    round(leer(descomposicion, "resolucion",    "Modelo"),  6),
    round(leer(descomposicion, "incertidumbre", "Modelo"),  6)
  )
) |>
  mutate(
    `Modelo − Mercado (pre-cierre)` =
      ifelse(is.na(Modelo), NA, round(Modelo - `Mercado (pre-cierre)`, 6))
  )

# Nota que viaja con el CSV: la fila REL no admite comparación entre
# competidores a este tamaño de muestra (ver Sección 2.1).
tabla_maestra <- tabla_maestra |>
  mutate(
    Nota = case_when(
      grepl("REL", Metrica) ~ paste(
        "NO comparable entre competidores: el Script 08 muestra que las tres",
        "REL están por debajo del percentil 95 del piso de ruido de binning."),
      grepl("UNC", Metrica) ~ "Depende solo de los datos observados, no del predictor.",
      TRUE                  ~ ""
    )
  )

cat("\n=== TABLA MAESTRA — COMPARACIÓN MODELO VS. MERCADO (24/25) ===\n")
print(tabla_maestra |> select(-Nota), n = Inf, width = Inf)

cat("\nNOTA sobre la fila REL: el Script 08 contrastó las tres calibraciones\n")
cat("contra un piso de ruido simulado y ninguna lo supera. La diferencia de\n")
cat("REL entre competidores NO es interpretable a este tamaño de muestra: la\n")
cat("brecha de Brier Score viene de RESOLUCIÓN. Ver la Sección 2.1.\n")

write.csv(tabla_maestra, "tabla_maestra_resultados.csv", row.names = FALSE)
cat("\nArchivo 'tabla_maestra_resultados.csv' guardado.\n")


# -----------------------------------------------------------------------------
# 2.1 CONTRASTES PRINCIPALES (Scripts 07 y 08)
# -----------------------------------------------------------------------------
# La tabla maestra reporta MAGNITUDES, no inferencia. Los contrastes que
# sostienen las afirmaciones del capítulo de resultados los produjeron los
# Scripts 07 y 08, y se exportan acá: sin ellos la síntesis entregaría
# diferencias puntuales sin ningún respaldo estadístico, que es lo primero que
# se pregunta al leer la tabla.
#
# Cada bloque se saltea si el objeto correspondiente no está en el .rds.

cat("\n", strrep("=", 74), "\n", sep = "")
cat("CONTRASTES PRINCIPALES — LOS NÚMEROS QUE SOSTIENEN LAS CONCLUSIONES\n")
cat(strrep("=", 74), "\n", sep = "")

# --- (a) Calibración: contraste primario del trabajo ---
if (!is.null(calibracion$contraste_primario)) {

  primario <- calibracion$contraste_primario

  cat("\n--- (a) Modelo vs. Mercado pre-cierre en Brier Score (Script 08) ---\n")
  cat(sprintf("Diferencia media (modelo − mercado): %+.6f\n", primario$dif_media))
  cat(sprintf("IC 95%%: [%+.6f ; %+.6f]\n", primario$ic95_inf, primario$ic95_sup))
  cat(sprintf("p = %.6f | Wilcoxon p = %.6f | d de Cohen = %.3f | n = %d\n",
              primario$p_t, primario$p_wilcoxon, primario$d_cohen,
              primario$n_pares))
  cat("Diferencia POSITIVA = el mercado predice mejor.\n")
  cat("Este es el contraste pre-especificado del trabajo y el que debe\n")
  cat("encabezar el capítulo de resultados.\n")

  write.csv(primario, "tabla_contraste_primario.csv", row.names = FALSE)
}

if (!is.null(calibracion$tests_significacion)) {
  write.csv(calibracion$tests_significacion,
            "tabla_contrastes_calibracion.csv", row.names = FALSE)
  cat("\nLos seis contrastes (3 comparaciones x 2 métricas) se exportan en\n")
  cat("'tabla_contrastes_calibracion.csv'.\n")
}

# --- (b) Sensibilidad a los equipos con parámetros imputados ---
if (!is.null(calibracion$contraste_primario_sin_asc)) {
  cat("\n--- (b) El mismo contraste excluyendo a los ascendidos (Script 08) ---\n")
  print(calibracion$contraste_primario_sin_asc |>
          select(n_pares, dif_media, ic95_inf, ic95_sup, p_t) |>
          mutate(across(where(is.numeric), ~ signif(., 4))))
  cat("La conclusión no depende de la imputación de parámetros.\n")
  write.csv(calibracion$contraste_primario_sin_asc,
            "tabla_contraste_sin_ascendidos.csv", row.names = FALSE)
}

# --- (c) Calibración: piso de ruido de REL ---
if (!is.null(calibracion$piso_ruido_rel)) {
  cat("\n--- (c) Calibración contra el piso de ruido de binning (Script 08) ---\n")
  print(calibracion$piso_ruido_rel |>
          mutate(across(where(is.numeric), ~ round(., 6))))
  if (!any(calibracion$piso_ruido_rel$excede_ruido)) {
    cat("Ninguno supera el percentil 95: los tres competidores son\n")
    cat("indistinguibles de un predictor perfectamente calibrado. La brecha de\n")
    cat("Brier Score es enteramente de RESOLUCIÓN, no de calibración.\n")
  }
  write.csv(calibracion$piso_ruido_rel, "tabla_piso_ruido_rel.csv",
            row.names = FALSE)
}

# --- (d) FLB: el test agregado, que es el instrumento de referencia ---
if (!is.null(mercado$flb_agregado)) {

  cat("\n--- (d) TEST AGREGADO DE FLB (Script 07, Sección 7.2) ---\n")
  cat("Interceptos por categoría + pendiente común | EE agrupados por partido\n")
  cat("H0: b = 1 (eficiencia) | FLB implica b > 1\n\n")

  print(mercado$flb_agregado |>
          select(serie, partidos, pendiente, se_agrupado,
                 ic95_inf, ic95_sup, z, p_valor) |>
          mutate(across(where(is.numeric), ~ round(., 4))),
        width = Inf)

  cat("\nESTE es el número que encabeza el capítulo de FLB, NO las seis\n")
  cat("regresiones desagregadas de la Sección 3. Aquella tabla incluye una\n")
  cat("fila que declara 'b > 1: consistente con FLB' —Visitante en la serie de\n")
  cat("cierre— y el Script 07 advierte expresamente que no corresponde\n")
  cat("encabezar con ese rechazo aislado: los seis contrastes comparten\n")
  cat("partidos y categorías, y el agregado es el que sintetiza la evidencia\n")
  cat("con inferencia válida.\n")

  write.csv(mercado$flb_agregado, "tabla_flb_agregado.csv", row.names = FALSE)
}

# --- (e) FLB en términos económicos: test de tendencia ---
if (!is.null(mercado$tendencia_cain)) {
  cat("\n--- (e) Test de tendencia del retorno según el precio (Script 07) ---\n")
  cat("retorno = alfa + beta * precio | EE agrupados por partido | H0: beta = 0\n\n")
  print(mercado$tendencia_cain |>
          mutate(across(where(is.numeric), ~ round(., 4))), width = Inf)
  cat("\nEs el contraste económico de referencia: usa las 4.560 observaciones\n")
  cat("individuales. El de bandas extremas se conserva por comparabilidad con\n")
  cat("Cain, Law y Peel (2000), pero descarta las bandas intermedias.\n")
  write.csv(mercado$tendencia_cain, "tabla_flb_tendencia.csv", row.names = FALSE)
}

if (!is.null(mercado$contraste_extremos)) {
  write.csv(mercado$contraste_extremos, "tabla_retornos_extremos.csv",
            row.names = FALSE)
}

# --- (f) Overround: robustez temporal del contraste pareado ---
if (!is.null(mercado$robustez_temporal)) {
  cat("\n--- (f) Overround pareado: robustez temporal (Script 07) ---\n")
  print(mercado$robustez_temporal |>
          mutate(across(where(is.numeric), ~ round(., 5))), width = Inf)
  cat("\nIMPRESCINDIBLE junto al contraste global de la Sección 3. El Script 07\n")
  cat("estableció que el efecto NO sobrevive al excluir 2021/22: el promedio de\n")
  cat("las cuatro temporadas mezcla dos regímenes. NO redactar 'el margen se\n")
  cat("expande hacia el cierre' como rasgo del período; el hallazgo es el\n")
  cat("CAMBIO DE RÉGIMEN.\n")
  write.csv(mercado$robustez_temporal, "tabla_overround_robustez_temporal.csv",
            row.names = FALSE)
}

if (!is.null(mercado$robustez_pareado)) {
  write.csv(mercado$robustez_pareado, "tabla_overround_robustez_extremos.csv",
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 3. TABLAS DEL ANÁLISIS DE MERCADO (CUATRO TEMPORADAS)
# -----------------------------------------------------------------------------
# Estas tablas corresponden al análisis del mercado, que NO involucra al modelo
# y por lo tanto usa la ventana completa de cuatro temporadas. Se exportan por
# separado, y el texto debe enunciar la regla de ventana temporal como criterio
# metodológico, no dejar que aparezca como una inconsistencia entre tablas.

# --- Overround por temporada y comparación pareada entre series ---
write.csv(mercado$overround_global,    "tabla_overround_global.csv",    row.names = FALSE)
write.csv(mercado$overround_temporada, "tabla_overround_temporada.csv", row.names = FALSE)
write.csv(mercado$pareado_temporada,   "tabla_overround_pareado.csv",   row.names = FALSE)

cat("\n=== OVERROUND — ANÁLISIS DE MERCADO (4 temporadas) ===\n")
print(mercado$overround_global)
cat("\nPor temporada:\n")
print(mercado$overround_temporada)

cat("\n--- Comparación pareada global (cierre − pre-cierre) ---\n")
cat(sprintf("  Diferencia media: %+.4f pp | t = %.4f | p = %.6f\n",
            mercado$test_t_pareado$estimate * 100,
            mercado$test_t_pareado$statistic,
            mercado$test_t_pareado$p.value))
cat(sprintf("  IC 95%%: [%+.4f ; %+.4f] pp | Wilcoxon p = %.6f\n",
            mercado$test_t_pareado$conf.int[1] * 100,
            mercado$test_t_pareado$conf.int[2] * 100,
            mercado$test_w_pareado$p.value))

# --- Regresiones de calibración: el resultado central del test de FLB ---
tabla_flb <- mercado$flb_regresiones |>
  select(serie, etiqueta_cat, n,
         intercepto, se_intercepto,
         pendiente, se_pendiente, ic95_inf, ic95_sup,
         z_pendiente, p_pendiente,
         wald_conjunto, p_conjunto, conclusion) |>
  mutate(across(where(is.numeric), ~ round(., 4)))

cat("\n=== REGRESIONES DE CALIBRACIÓN — TEST DE FLB ===\n")
print(tabla_flb, n = Inf, width = Inf)

write.csv(tabla_flb, "tabla_flb_regresiones.csv", row.names = FALSE)

# --- Potencia del test ---
tabla_potencia_exp <- mercado$tabla_potencia |>
  select(etiqueta_cat, n, sd_logit, b_alt, potencia) |>
  mutate(across(where(is.numeric), ~ round(., 4)))

cat("\n=== POTENCIA DEL TEST DE FLB ===\n")
print(tabla_potencia_exp, n = Inf, width = Inf)

write.csv(tabla_potencia_exp, "tabla_flb_potencia.csv", row.names = FALSE)

# --- Retornos por banda de precio (Cain, Law y Peel, 2000) ---
write.csv(mercado$tabla_cain, "tabla_retornos_por_banda.csv", row.names = FALSE)

cat("\n=== RETORNOS POR BANDA DE PRECIO (Cain, Law y Peel, 2000) ===\n")
print(mercado$tabla_cain |> mutate(across(where(is.numeric), ~ round(., 4))),
      n = Inf, width = Inf)

# --- Correlación descriptiva (continuidad con la versión anterior) ---
write.csv(mercado$correlaciones, "tabla_correlacion_descriptiva.csv", row.names = FALSE)

cat("\nTablas del análisis de mercado exportadas.\n")


# -----------------------------------------------------------------------------
# 4. GRÁFICOS DEL SCRIPT 07 — EXPORTADOS COMO .PNG
# -----------------------------------------------------------------------------
# Los gráficos del análisis de mercado se regeneran acá con formato de
# exportación. Se levantan los objetos del Script 07 en lugar de recalcularlos,
# para que no puedan divergir de lo que se reportó allá.

datos_mercado_4t <- mercado$datos_mercado
mercado_largo    <- mercado$mercado_largo
bins_mercado     <- mercado$bins_mercado

col_series <- c("Pre-cierre" = "#2c3e50", "Cierre" = "#16a085")

# --- Gráfico 07-A: distribución del overround, ambas series ---
g07a <- datos_mercado_4t |>
  select(overround, overround_c) |>
  pivot_longer(everything(), names_to = "serie", values_to = "valor") |>
  filter(!is.na(valor)) |>
  mutate(serie = if_else(serie == "overround",
                         "Pre-cierre (B365H/D/A)", "Cierre (B365CH/CD/CA)")) |>
  ggplot(aes(x = valor * 100, fill = serie)) +
  geom_histogram(binwidth = 0.2, color = "white", alpha = 0.75,
                 position = "identity") +
  facet_wrap(~serie, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Pre-cierre (B365H/D/A)" = "#2c3e50",
                               "Cierre (B365CH/CD/CA)"  = "#16a085")) +
  labs(
    title    = "Distribución del overround — Bet365",
    subtitle = sprintf("Premier League 2021/22–2024/25 | N = %d partidos",
                       nrow(datos_mercado_4t)),
    x = "Overround (%)", y = "Frecuencia", fill = NULL
  ) +
  theme_minimal() + theme(legend.position = "none")

ggsave("grafico_07a_overround_distribucion.png", g07a, width = 8, height = 6, dpi = 300)

# --- Gráfico 07-B: diagrama de fiabilidad con IC de Wilson ---
g07b <- bins_mercado |>
  ggplot(aes(x = pi_med, y = frec_real, color = serie)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbar(aes(ymin = inf, ymax = sup), width = 0.012, alpha = 0.7) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  scale_color_manual(values = col_series) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Diagrama de fiabilidad del mercado",
    subtitle = "Deciles con IC de Wilson al 95% | 2021/22–2024/25 | Diagonal = calibración perfecta",
    x = "Probabilidad implícita normalizada (promedio del decil)",
    y = "Frecuencia real de ocurrencia", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_07b_fiabilidad_deciles.png", g07b, width = 8, height = 5, dpi = 300)

# --- Gráfico 07-C: pendientes de la regresión de calibración con IC ---
datos_07c <- mercado$flb_regresiones |>
  mutate(
    # Invertido: con el giro, el primer nivel queda abajo. Así "Local (H)"
    # aparece arriba y se conserva el orden H, D, A del resto del trabajo.
    etiqueta_cat = factor(etiqueta_cat,
                          levels = c("Visitante (A)", "Empate (D)", "Local (H)")),
    # La serie primaria va primero en la leyenda. Si la etiqueta de la serie
    # no fuera exactamente "Pre-cierre", ajustar el literal.
    serie = factor(serie, levels = c("Pre-cierre", "Cierre"))
  )

g07c <- ggplot(datos_07c,
               aes(x = pendiente, y = etiqueta_cat, color = serie)) +
  
  # Referencia de eficiencia.
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray45", linewidth = 0.6) +
  
  # Magnitud que los intervalos descartan. Es el resultado que el texto
  # reporta en formulación acotada; conviene que la figura lo exhiba.
  geom_vline(xintercept = 1.30, linetype = "dotted",
             color = "gray65", linewidth = 0.6) +
  
  geom_pointrange(aes(xmin = ic95_inf, xmax = ic95_sup),
                  position = position_dodge(width = 0.55),
                  size = 0.55, linewidth = 0.7) +
  
  scale_x_continuous(
    breaks = c(0.6, 0.8, 1.0, 1.2, 1.4),
    labels = c("0,6", "0,8", "1,0", "1,2", "1,4"),
    sec.axis = dup_axis(breaks = c(1, 1.30),
                        labels = c("eficiencia", "1,3"),
                        name = NULL)) +
  
  labs(
    title    = "Pendiente de la regresión de calibración, por categoría y serie",
    subtitle = paste0("Intervalos de confianza del 95 % | la eficiencia implica ",
                      "b = 1 y el sesgo b > 1"),
    x = "Pendiente estimada (b)", y = NULL, color = NULL
  ) +
  
  theme_minimal(base_size = 10, base_family = "serif") +
  theme(
    legend.position   = "bottom",
    legend.margin     = margin(t = -4),
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title        = element_text(size = 12),
    plot.subtitle     = element_text(color = "gray35", size = 8.5,
                                     margin = margin(b = 8)),
    axis.title.x      = element_text(size = 9),
    axis.text.y       = element_text(size = 10),
    axis.text.x.top   = element_text(color = "gray45", size = 8)
  )

ggsave("grafico_07c_flb_pendientes.png", g07c,
       width = 8, height = 3.4, dpi = 300)

# --- Gráfico 07-D: overround según probabilidad del favorito ---
g07d <- datos_mercado_4t |>
  mutate(pi_max = pmax(pi_H, pi_A), bin_mercado = ntile(pi_max, 10)) |>
  summarise(
    pi_max_med      = mean(pi_max),
    overround_med   = mean(overround,   na.rm = TRUE) * 100,
    overround_c_med = mean(overround_c, na.rm = TRUE) * 100,
    .by = bin_mercado
  ) |>
  pivot_longer(cols = c(overround_med, overround_c_med),
               names_to = "serie", values_to = "overround") |>
  mutate(serie = if_else(serie == "overround_med", "Pre-cierre", "Cierre")) |>
  ggplot(aes(x = pi_max_med * 100, y = overround, color = serie)) +
  geom_point(size = 3.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_color_manual(values = col_series) +
  labs(
    title    = "Overround según la probabilidad del favorito",
    subtitle = "Cada punto es un decil de partidos | Premier League 2021/22–2024/25",
    x = "Probabilidad implícita del favorito (promedio del decil, %)",
    y = "Overround promedio (%)", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_07d_overround_vs_favorito.png", g07d, width = 8, height = 5, dpi = 300)

# --- Gráfico 07-E: evolución del overround por temporada ---
g07e <- mercado$overround_temporada |>
  select(temporada, pre_cierre, cierre) |>
  pivot_longer(cols = c(pre_cierre, cierre),
               names_to = "serie", values_to = "overround") |>
  mutate(serie = if_else(serie == "pre_cierre", "Pre-cierre", "Cierre")) |>
  ggplot(aes(x = temporada, y = overround, color = serie, group = serie)) +
  geom_line(linewidth = 0.9) + geom_point(size = 3.5) +
  scale_color_manual(values = col_series) +
  labs(
    title    = "Evolución del overround medio por temporada",
    subtitle = "Bet365 — Premier League",
    x = NULL, y = "Overround medio (%)", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_07e_overround_temporada.png", g07e, width = 8, height = 5, dpi = 300)

# --- Gráfico 07-F: retornos por banda de precio (Cain, Law y Peel, 2000) ---
# LA LÍNEA DE REFERENCIA NO ES CERO. Con overround ninguna apuesta puede rendir
# cero: bajo calibración perfecta y margen uniforme, CUALQUIER apuesta unitaria
# rinde 1/(1+OR) − 1, alrededor de −5,1%. Una línea en cero haría parecer que
# las cuatro bandas pierden dinero, cuando lo que hay que juzgar es su desvío
# respecto de esa referencia. Es la corrección (g) del Script 07, que esta
# versión de exportación deshacía.
datos_g07f <- if (!is.null(mercado$tabla_cain_ref)) {
  mercado$tabla_cain_ref
} else {
  mercado$tabla_cain |>
    left_join(mercado$referencia_justa |> select(serie, retorno_justo),
              by = "serie")
}

g07f <- datos_g07f |>
  ggplot(aes(x = banda, y = retorno_medio, color = serie, group = serie)) +
  geom_hline(aes(yintercept = retorno_justo, color = serie),
             linetype = "dashed", linewidth = 0.7, alpha = 0.8) +
  geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup),
                  position = position_dodge(width = 0.3), size = 0.6) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.6, alpha = 0.6) +
  scale_color_manual(values = col_series) +
  labs(
    title    = "Retorno a una apuesta unitaria por banda de precio",
    subtitle = paste0(
      "Réplica de Cain, Law y Peel (2000, Tabla 1) | IC del 95% | 2021/22-2024/25\n",
      "Línea discontinua: retorno esperado bajo calibración perfecta y margen ",
      "uniforme, 1/(1+OR) - 1"),
    x = "Banda de precio i = 1/cuota",
    y = "Retorno medio por unidad apostada", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_07f_retornos_por_banda.png", g07f, width = 8, height = 5, dpi = 300)

# --- Gráfico 07-G: curva de potencia del test de FLB ---
g07g <- mercado$tabla_potencia |>
  mutate(etiqueta_cat = factor(etiqueta_cat,
                               levels = c("Local (H)", "Empate (D)", "Visitante (A)"))) |>
  ggplot(aes(x = b_alt, y = potencia, color = etiqueta_cat, group = etiqueta_cat)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.9) + geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = c("Local (H)" = "#2980b9",
                                "Empate (D)" = "#e67e22",
                                "Visitante (A)" = "#27ae60")) +
  labs(
    title    = "Potencia del test de FLB por categoría",
    subtitle = "Simulación sobre la distribución observada de π | α = 0,05 | Línea: potencia del 80%",
    x = "Tamaño del sesgo bajo la alternativa (pendiente verdadera b)",
    y = "Potencia", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_07g_potencia_flb.png", g07g, width = 8, height = 5, dpi = 300)

cat("Gráficos del Script 07 exportados (07a a 07g).\n")

# -----------------------------------------------------------------------------
# 5. GRÁFICOS DEL SCRIPT 08 — EXPORTADOS COMO .PNG
# -----------------------------------------------------------------------------

col_comp <- c("Mercado (pre-cierre)" = "#e74c3c",
              "Mercado (cierre)"     = "#16a085",
              "Modelo Maher"         = "#2980b9")

renombrar_bs <- function(x) case_when(x == "bs_mercado" ~ "Mercado (pre-cierre)",
                                      x == "bs_cierre"  ~ "Mercado (cierre)",
                                      TRUE              ~ "Modelo Maher")
renombrar_ll <- function(x) case_when(x == "ll_mercado" ~ "Mercado (pre-cierre)",
                                      x == "ll_cierre"  ~ "Mercado (cierre)",
                                      TRUE              ~ "Modelo Maher")

# --- Gráfico 08-A: Brier Score por decil ---
g08a <- metricas_por_bin |>
  select(bin, bs_mercado, bs_cierre, bs_modelo) |>
  pivot_longer(cols = c(bs_mercado, bs_cierre, bs_modelo),
               names_to = "fuente", values_to = "brier") |>
  mutate(fuente = renombrar_bs(fuente)) |>
  ggplot(aes(x = factor(bin), y = brier, color = fuente, group = fuente)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  scale_color_manual(values = col_comp) +
  labs(
    title    = "Brier Score por decil de probabilidad",
    subtitle = "Menor es mejor | Deciles definidos sobre la serie pre-cierre | Premier League 2024/25",
    x = "Decil de probabilidad implícita (1 = longshots, 10 = favoritos)",
    y = "Brier Score", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_08a_brier_deciles.png", g08a,
       width = 8, height = 5, dpi = 300)
cat("Guardado: grafico_08a_brier_deciles.png\n")


# --- Gráfico 08-B: log-loss por decil ---
g08b <- metricas_por_bin |>
  select(bin, ll_mercado, ll_cierre, ll_modelo) |>
  pivot_longer(cols = c(ll_mercado, ll_cierre, ll_modelo),
               names_to = "fuente", values_to = "logloss") |>
  mutate(fuente = renombrar_ll(fuente)) |>
  ggplot(aes(x = factor(bin), y = logloss, color = fuente, group = fuente)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  scale_color_manual(values = col_comp) +
  labs(
    title    = "Log-loss condicional por decil de probabilidad",
    subtitle = paste0(
      "Premier League 2024/25 | Calculado solo sobre los resultados que ocurrieron\n",
      "DESCRIPTIVO: al condicionar sobre el resultado deja de ser una regla de ",
      "puntuación propia y NO ordena competidores"),
    x = "Decil de probabilidad implícita (1 = longshots, 10 = favoritos)",
    y = "Log-loss condicional", color = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_08b_logloss_deciles.png", g08b,
       width = 8, height = 5, dpi = 300)
cat("Guardado: grafico_08b_logloss_deciles.png\n")


# --- Gráfico 08-C: descomposición de Murphy ---

# Formateador con coma decimal, para etiquetas y eje.
fmt_coma <- function(x, dig = 4) {
  formatC(x, format = "f", digits = dig, decimal.mark = ",")
}

g08c <- tibble(
  componente = factor(rep(c("Calibración (REL)", "Resolución (RES)"), 3),
                      levels = c("Calibración (REL)", "Resolución (RES)")),
  fuente = factor(c(rep("Mercado (pre-cierre)", 2),
                    rep("Mercado (de cierre)",  2),
                    rep("Modelo Maher",         2)),
                  levels = c("Mercado (pre-cierre)",
                             "Mercado (de cierre)",
                             "Modelo Maher")),
  valor  = c(
    leer(descomposicion, "calibracion", "Mercado (pre-cierre)"),
    leer(descomposicion, "resolucion",  "Mercado (pre-cierre)"),
    leer(descomposicion, "calibracion", "Mercado (cierre)"),
    leer(descomposicion, "resolucion",  "Mercado (cierre)"),
    leer(descomposicion, "calibracion", "Modelo"),
    leer(descomposicion, "resolucion",  "Modelo")
  )
) |>
  ggplot(aes(x = componente, y = valor, fill = fuente)) +
  
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9, width = 0.8) +
  
  geom_text(aes(label = fmt_coma(valor)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 2.6, family = "serif", color = "gray20") +
  
  # Los valores de calibración son pequeños frente a los de resolución; el
  # margen superior evita que las etiquetas queden pegadas al borde.
  scale_y_continuous(labels = function(x) fmt_coma(x, 2),
                     expand = expansion(mult = c(0, 0.12))) +
  
  # Los nombres de las series se fijan con `values` para no depender del orden
  # de los niveles; si col_comp usara las etiquetas anteriores, ajustarlo acá.
  scale_fill_manual(values = unname(col_comp)) +
  
  labs(
    title    = "Descomposición del Brier Score",
    subtitle = sprintf(
      paste0("En resolución, un valor mayor indica más capacidad de ",
             "discriminación entre partidos\n",
             "Las diferencias de calibración no son interpretables: las tres ",
             "se ubican dentro del piso de ruido\n",
             "La incertidumbre, %s, es idéntica para los tres y se omite"),
      fmt_coma(leer(descomposicion, "incertidumbre", "Modelo"))),
    x = NULL, y = "Valor", fill = NULL
  ) +
  
  theme_minimal(base_size = 10, base_family = "serif") +
  theme(
    legend.position   = "bottom",
    legend.margin     = margin(t = -4),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title        = element_text(size = 12),
    plot.subtitle     = element_text(color = "gray35", size = 8.5,
                                     margin = margin(b = 8)),
    axis.title.y      = element_text(size = 9),
    axis.text.x       = element_text(size = 10)
  )

ggsave("grafico_08c_descomposicion_murphy.png", g08c,
       width = 7, height = 3.6, dpi = 300)
cat("Guardado: grafico_08c_descomposicion_murphy.png\n")


# --- Gráfico 08-D: Brier Score por tipo de resultado ---
g08d <- metricas_por_resultado |>
  select(resultado, bs_mercado, bs_cierre, bs_modelo) |>
  pivot_longer(cols = c(bs_mercado, bs_cierre, bs_modelo),
               names_to = "fuente", values_to = "brier") |>
  mutate(
    fuente    = renombrar_bs(fuente),
    resultado = factor(resultado, levels = c("H", "D", "A"),
                       labels = c("Local (H)", "Empate (D)", "Visitante (A)"))
  ) |>
  ggplot(aes(x = resultado, y = brier, fill = fuente)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = col_comp) +
  labs(
    title    = "Brier Score por tipo de resultado",
    subtitle = "Menor es mejor | Premier League 2024/25",
    x = NULL, y = "Brier Score", fill = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("grafico_08d_brier_por_resultado.png", g08d,
       width = 8, height = 5, dpi = 300)
cat("Guardado: grafico_08d_brier_por_resultado.png\n")


# -----------------------------------------------------------------------------
# 6. GRÁFICO DE SÍNTESIS: MODELO VS. MERCADO PARTIDO A PARTIDO
# -----------------------------------------------------------------------------
# Cada punto es un partido. Eje X = probabilidad del mercado para el resultado
# que efectivamente ocurrió; eje Y = probabilidad del modelo para ese mismo
# resultado. La diagonal es la línea de acuerdo perfecto entre ambos.
# Puntos por encima: el modelo fue más optimista que el mercado para ese resultado.
# Puntos por debajo: el modelo fue más pesimista.
# Coloreamos por tipo de resultado (H, D, A).

datos_resultado_ocurrido <- datos_largo |>
  filter(ocurrio == 1)   # solo la fila del resultado que realmente ocurrió

g_sintesis <- datos_resultado_ocurrido |>
  mutate(
    resultado = factor(resultado,
                       levels = c("H", "D", "A"),
                       labels = c("Local (H)", "Empate (D)", "Visitante (A)"))
  ) |>
  ggplot(aes(x = pi_mercado, y = prob_modelo, color = resultado)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_point(alpha = 0.5, size = 1.8) +
  scale_color_manual(values = c(
    "Local (H)"    = "#2980b9",
    "Empate (D)"   = "#e67e22",
    "Visitante (A)"= "#27ae60"
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title    = "Modelo vs. mercado: probabilidad del resultado ocurrido",
    subtitle = "Cada punto es un partido | Diagonal = acuerdo perfecto | Premier League 2024/25",
    x        = "Probabilidad implícita del mercado (Bet365)",
    y        = "Probabilidad del modelo desarrollado",
    color    = "Resultado"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("grafico_09a_sintesis_modelo_vs_mercado.png", g_sintesis,
       width = 8, height = 7, dpi = 300)
cat("Guardado: grafico_09a_sintesis_modelo_vs_mercado.png\n")


# -----------------------------------------------------------------------------
# 7. PARTIDOS CON MAYOR DISCREPANCIA MODELO VS. MERCADO
# -----------------------------------------------------------------------------
# Para cada partido calculamos la discrepancia total como la suma de las
# diferencias absolutas entre probabilidades del modelo y del mercado
# para los tres resultados posibles.
# Los 10 partidos con mayor discrepancia son los casos donde el modelo
# y el mercado disentían más en su evaluación del partido.

discrepancia <- datos |>
  mutate(
    disc_H     = abs(prob_H - pi_H),
    disc_D     = abs(prob_D - pi_D),
    disc_A     = abs(prob_A - pi_A),
    disc_total = disc_H + disc_D + disc_A
  ) |>
  arrange(desc(disc_total)) |>
  select(
    date, home_team, away_team, ftr,
    prob_H, pi_H,
    prob_D, pi_D,
    prob_A, pi_A,
    disc_total
  ) |>
  mutate(across(where(is.numeric), ~ round(., 4)))

top10_mayor <- slice_head(discrepancia, n = 10)
top10_menor <- slice_tail(discrepancia, n = 10)

cat("\n=== TOP 10 PARTIDOS: MAYOR DISCREPANCIA MODELO VS. MERCADO ===\n")
print(top10_mayor, n = 10)

cat("\n=== TOP 10 PARTIDOS: MENOR DISCREPANCIA MODELO VS. MERCADO ===\n")
print(top10_menor, n = 10)

write.csv(top10_mayor, "tabla_discrepancia_top10_mayor.csv", row.names = FALSE)
write.csv(top10_menor, "tabla_discrepancia_top10_menor.csv", row.names = FALSE)
cat("\nArchivos de discrepancia guardados.\n")


# -----------------------------------------------------------------------------
# 8. PARÁMETROS DEL MODELO FINAL
# -----------------------------------------------------------------------------
# Los parámetros de ataque y defensa estimados para cada equipo son un
# resultado sustantivo de la tesina: permiten ver qué equipos el modelo
# identifica como más fuertes ofensiva y defensivamente.
# Levantamos los parámetros desde resultados_modelo_final.rds.

# VENTANA TEMPORAL DE LOS PARÁMETROS. Se estiman sobre 2021/22-2023/24 y se
# APLICAN a la validación 2024/25. Por eso la tabla incluye equipos que no
# disputaron la 24/25: descendieron en alguna de las temporadas de
# entrenamiento. Rotularla como "Premier League 2024/25" era incorrecto.
#
# Y NO TODOS LOS PARÁMETROS ESTÁN ESTIMADOS. Los equipos ascendidos a 24/25 no
# tienen historial en la muestra de entrenamiento: sus valores son IMPUTADOS
# vía peso_descendidos (Script 04). Mostrarlos junto a los estimados sin
# distinguirlos induce a error, sobre todo después del análisis de sensibilidad
# de la Sección 4.6 del Script 08.

ataques  <- modelo_final$ataques_fin
defensas <- modelo_final$defensas_fin

# Equipos que efectivamente disputaron la temporada de validación
equipos_2425 <- sort(unique(c(datos$home_team, datos$away_team)))

# Equipos con parámetros imputados (los guarda el Script 08)
equipos_imputados <- if (!is.null(calibracion$equipos_nuevos)) {
  calibracion$equipos_nuevos
} else {
  cat("AVISO: no está 'equipos_nuevos' en resultados_calibracion.rds;\n")
  cat("       no se pueden marcar los parámetros imputados.\n")
  character(0)
}

tabla_parametros <- tibble(
  equipo  = names(ataques),
  ataque  = round(as.numeric(ataques),  4),
  defensa = round(as.numeric(defensas), 4)
) |>
  mutate(
    imputado  = equipo %in% equipos_imputados,
    jugo_2425 = equipo %in% equipos_2425
  ) |>
  arrange(desc(ataque))

cat("\n=== PARÁMETROS DEL MODELO FINAL (ordenados por ataque) ===\n")
cat("Estimados sobre 2021/22-2023/24. 'imputado' = sin historial de\n")
cat("entrenamiento; 'jugo_2425' = disputó la temporada de validación.\n\n")
print(tabla_parametros, n = Inf)

cat(sprintf("\nEquipos en la tabla: %d | disputaron 24/25: %d | imputados: %d\n",
            nrow(tabla_parametros), sum(tabla_parametros$jugo_2425),
            sum(tabla_parametros$imputado)))

if (any(!tabla_parametros$jugo_2425)) {
  cat("No disputaron la 24/25 (descendidos en temporadas de entrenamiento): ",
      paste(tabla_parametros$equipo[!tabla_parametros$jugo_2425],
            collapse = ", "), "\n", sep = "")
  cat("Sus parámetros son un resultado legítimo de la estimación, pero NO\n")
  cat("intervienen en la validación 24/25. Declararlo si la tabla va al texto.\n")
}

if (any(tabla_parametros$imputado)) {
  cat("Con parámetros IMPUTADOS (no estimados): ",
      paste(tabla_parametros$equipo[tabla_parametros$imputado],
            collapse = ", "), "\n", sep = "")
}

write.csv(tabla_parametros, "tabla_parametros_modelo.csv", row.names = FALSE)
cat("Archivo 'tabla_parametros_modelo.csv' guardado.\n")

library(ggrepel)

# --- Interruptor de proporción ------------------------------------------------
# TRUE  -> una unidad de ataque mide lo mismo que una de defensa. La diagonal
#          queda a 45 grados y comparar ambas dimensiones es legítimo. El panel
#          toma el aspecto que le dan los datos y no se lo puede ensanchar.
# FALSE -> el panel se estira al tamaño del dispositivo y admite un formato
#          apaisado. La diagonal se retira sola, porque fuera de proporción
#          sugeriría magnitudes que no representa.

proporcion_fija <- TRUE

# --- Preparación de los datos -------------------------------------------------
datos_grafico <- tabla_parametros |>
  mutate(
    caso = case_when(
      imputado   ~ "Parámetros imputados",
      !jugo_2425 ~ "No disputó 2024/25",
      TRUE       ~ "Estimado y evaluado"
    ),
    caso = factor(caso, levels = c("Estimado y evaluado",
                                   "No disputó 2024/25",
                                   "Parámetros imputados"))
  )

# --- Límites ------------------------------------------------------------------
# Se toman del rango efectivo de los datos y no de un máximo simétrico. Esto es
# lo que elimina el espacio vertical sobrante de la versión anterior, donde los
# límites quedaban gobernados por el ataque de Manchester City.

margen  <- 0.06
rango_x <- range(datos_grafico$ataque)  + c(-margen, margen)
rango_y <- range(datos_grafico$defensa) + c(-margen, margen)

# --- Gráfico ------------------------------------------------------------------
g_parametros <- ggplot(datos_grafico, aes(x = ataque, y = defensa)) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray55", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray55", linewidth = 0.4) +
  
  geom_point(aes(shape = caso, fill = caso),
             size = 2.6, color = "gray25", stroke = 0.4, alpha = 0.9) +
  
  # Etiquetas de los 26 equipos. ggrepel las desplaza hasta que no se solapan y
  # traza un segmento fino hacia el punto solo cuando debe alejarlas.
  geom_text_repel(
    aes(label = equipo),
    family             = "serif",
    size               = 3.0,
    color              = "gray20",
    segment.color      = "gray65",
    segment.size       = 0.25,
    min.segment.length = 0.25,   # traza segmento recién a partir de esta distancia
    box.padding        = 0.40,   # espacio libre alrededor de cada etiqueta
    point.padding      = 0.20,   # espacio libre alrededor de cada punto
    max.overlaps       = Inf,    # nunca descartar una etiqueta
    force              = 1.5,    # intensidad de la repulsión
    max.iter           = 20000,
    seed               = 2024    # fija el resultado: el gráfico es reproducible
  ) +
  
  scale_shape_manual(values = c("Estimado y evaluado"  = 21,
                                "No disputó 2024/25"   = 24,
                                "Parámetros imputados" = 23)) +
  scale_fill_manual(values  = c("Estimado y evaluado"  = "#2c3e50",
                                "No disputó 2024/25"   = "#bdc3c7",
                                "Parámetros imputados" = "#f39c12")) +
  
  # El eje vertical se invierte: los valores negativos (mejor defensa) arriba.
  scale_y_reverse(limits = rev(rango_y),
                  expand = expansion(mult = 0.01),
                  labels = scales::label_number(decimal.mark = ",")) +
  scale_x_continuous(limits = rango_x,
                     expand = expansion(mult = 0.01),
                     labels = scales::label_number(decimal.mark = ",")) +
  
  labs(
    title    = "Ataque y defensa estimados, por equipo",
    subtitle = paste0(
      "Modelo de Maher estimado sobre 2021/22-2023/24, en escala logarítmica\n",
      "El origen es el equipo promedio de la liga; arriba y a la derecha, ",
      "mejor en ambas dimensiones"),
    x        = "Parámetro de ataque  (derecha = convierte más)",
    y        = "Parámetro de defensa  (arriba = recibe menos)",
    shape    = NULL,
    fill     = NULL
  ) +
  
  theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    legend.position  = "bottom",
    legend.margin    = margin(t = -4),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(color = "gray35", size = 8.5,
                                    margin = margin(b = 8)),
    axis.title       = element_text(size = 9)
  )

# La diagonal de equivalencia solo se dibuja con proporción fija.
if (proporcion_fija) {
  g_parametros <- g_parametros +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "gray78", linewidth = 0.4) +
    coord_fixed()
}

# --- Guardado -----------------------------------------------------------------
if (proporcion_fija) {
  ancho <- 7.0; alto <- 6.0
} else {
  ancho <- 8.5; alto <- 5.0
}

ggsave("grafico_09b_parametros_modelo.png", g_parametros,
       width = 6.5, height = 6.4, dpi = 300)
cat("Guardado: grafico_09b_parametros_modelo.png\n")


# -----------------------------------------------------------------------------
# 9. PARTIDOS DONDE EL MODELO ACERTÓ MÁS Y FALLÓ MÁS
# -----------------------------------------------------------------------------
# Para cada partido tomamos la probabilidad que el modelo asignó al resultado
# que realmente ocurrió. Cuanto más alta esa probabilidad, mejor fue la
# predicción del modelo para ese partido. Cuanto más baja, peor.
# Esto es distinto a la discrepancia con el mercado: acá el criterio es
# la calidad de la predicción del modelo contra la realidad observada.

aciertos_modelo <- datos_largo |>
  filter(ocurrio == 1) |>
  select(date, home_team, away_team, ftr,
         pi_mercado, prob_modelo) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  arrange(desc(prob_modelo))

top10_aciertos <- slice_head(aciertos_modelo, n = 10)
top10_fallos   <- slice_tail(aciertos_modelo, n = 10) |> arrange(prob_modelo)

cat("\n=== TOP 10 PARTIDOS: MODELO ASIGNÓ MAYOR PROBABILIDAD AL RESULTADO CORRECTO ===\n")
print(top10_aciertos, n = 10)

cat("\n=== TOP 10 PARTIDOS: MODELO ASIGNÓ MENOR PROBABILIDAD AL RESULTADO CORRECTO ===\n")
print(top10_fallos, n = 10)

write.csv(top10_aciertos, "tabla_modelo_top10_aciertos.csv", row.names = FALSE)
write.csv(top10_fallos,   "tabla_modelo_top10_fallos.csv",   row.names = FALSE)
cat("\nArchivos de aciertos y fallos del modelo guardados.\n")


# -----------------------------------------------------------------------------
# 10. ESTADÍSTICAS DESCRIPTIVAS DE LA TEMPORADA 24/25
# -----------------------------------------------------------------------------
# Contexto para el lector de la tesina: cómo fue la temporada de validación
# en términos de resultados, goles y equilibrio competitivo.

cat("\n=== ESTADÍSTICAS DESCRIPTIVAS — PREMIER LEAGUE 2024/25 ===\n")

# Distribución de resultados
dist_resultados <- datos |>
  count(ftr) |>
  mutate(
    resultado   = case_when(
      ftr == "H" ~ "Victoria local",
      ftr == "D" ~ "Empate",
      ftr == "A" ~ "Victoria visitante"
    ),
    porcentaje = round(n / sum(n) * 100, 1)
  ) |>
  select(resultado, n, porcentaje)

cat("\nDistribución de resultados:\n")
print(dist_resultados)

# Estadísticas de goles
stats_goles <- datos |>
  summarise(
    goles_totales       = sum(fthg + ftag),
    goles_por_partido   = round(mean(fthg + ftag), 3),
    goles_local_prom    = round(mean(fthg), 3),
    goles_visitante_prom = round(mean(ftag), 3),
    partidos_sin_goles  = sum(fthg + ftag == 0),
    partidos_mas_5_goles = sum(fthg + ftag > 5)
  )

cat("\nEstadísticas de goles:\n")
print(stats_goles)

# Overround: recordatorio del resultado clave
cat("\nOverround Bet365:\n")
datos |>
  summarise(
    media   = round(mean(overround) * 100, 3),
    mediana = round(median(overround) * 100, 3),
    desvio  = round(sd(overround) * 100, 3)
  ) |>
  print()

# ---------------------------------------------------------------------------
# % de partidos donde ganó el favorito del mercado
# Favorito = resultado con mayor probabilidad implícita normalizada
# Conecta con el análisis de calibración: si el mercado está bien calibrado,
# el favorito debería ganar en proporción similar a su probabilidad promedio.
# ---------------------------------------------------------------------------
pct_favorito_gana <- datos |>
  mutate(
    favorito = case_when(
      pi_H >= pi_D & pi_H >= pi_A ~ "H",
      pi_D >= pi_H & pi_D >= pi_A ~ "D",
      TRUE                         ~ "A"
    ),
    favorito_gano = as.integer(favorito == ftr)
  ) |>
  summarise(
    partidos             = n(),
    favorito_gano_n      = sum(favorito_gano),
    favorito_gano_pct    = round(mean(favorito_gano) * 100, 1),
    # Probabilidad media que el mercado asignó al favorito: es la referencia
    # contra la que hay que leer el porcentaje observado. Sin ella el bloque
    # anunciaba una comparación que no hacía.
    prob_media_favorito  = round(mean(pmax(pi_H, pi_D, pi_A)) * 100, 1)
  )

cat("\nResultados del favorito según mercado:\n")
print(pct_favorito_gana)
cat(sprintf("El favorito ganó el %.1f%% de los partidos; el mercado le asignaba\n",
            pct_favorito_gana$favorito_gano_pct))
cat(sprintf("en promedio un %.1f%% de probabilidad. Diferencia: %+.1f pp.\n",
            pct_favorito_gana$prob_media_favorito,
            pct_favorito_gana$favorito_gano_pct -
              pct_favorito_gana$prob_media_favorito))
cat("Es un descriptivo agregado, no un test de calibración: la inferencia\n")
cat("está en el Script 08 y en la Sección 2.1.\n")

# ---------------------------------------------------------------------------
# Goles promedio por mes (evolución temporal de la temporada)
# Permite detectar períodos atípicos que puedan haber afectado al modelo.
# ---------------------------------------------------------------------------
goles_por_mes <- datos |>
  mutate(mes = format(as.Date(date), "%Y-%m")) |>
  group_by(mes) |>
  summarise(
    partidos          = n(),
    goles_por_partido = round(mean(fthg + ftag), 3),
    .groups = "drop"
  ) |>
  arrange(mes)

cat("\nGoles promedio por mes:\n")
print(goles_por_mes, n = Inf)

# ---------------------------------------------------------------------------
# Distribución de marcadores finales más frecuentes (top 10)
# Conecta con Maher (1982): el modelo Poisson debería replicar esta distribución.
# Marcadores bajos (0-0, 1-0, 0-1, 1-1) son los que Dixon & Coles (1997)
# identificaron como sub/sobrerepresentados respecto al modelo independiente.
# ---------------------------------------------------------------------------
top_marcadores <- datos |>
  mutate(marcador = paste0(fthg, "-", ftag)) |>
  count(marcador, name = "frecuencia") |>
  arrange(desc(frecuencia)) |>
  slice_head(n = 10) |>
  mutate(porcentaje = round(frecuencia / nrow(datos) * 100, 1))

cat("\nTop 10 marcadores finales más frecuentes:\n")
print(top_marcadores)

# Exportar tabla descriptiva completa
tabla_descriptiva <- bind_rows(
  tibble(estadistica = "Partidos totales",                valor = as.character(nrow(datos))),
  tibble(estadistica = "Victorias locales (%)",           valor = as.character(dist_resultados$porcentaje[dist_resultados$resultado == "Victoria local"])),
  tibble(estadistica = "Empates (%)",                     valor = as.character(dist_resultados$porcentaje[dist_resultados$resultado == "Empate"])),
  tibble(estadistica = "Victorias visitante (%)",         valor = as.character(dist_resultados$porcentaje[dist_resultados$resultado == "Victoria visitante"])),
  tibble(estadistica = "Favorito del mercado ganó (%)",   valor = as.character(pct_favorito_gana$favorito_gano_pct)),
  tibble(estadistica = "Goles por partido",               valor = as.character(stats_goles$goles_por_partido)),
  tibble(estadistica = "Goles local promedio",            valor = as.character(stats_goles$goles_local_prom)),
  tibble(estadistica = "Goles visitante promedio",        valor = as.character(stats_goles$goles_visitante_prom)),
  tibble(estadistica = "Partidos sin goles",              valor = as.character(stats_goles$partidos_sin_goles)),
  tibble(estadistica = "Partidos con más de 5 goles",     valor = as.character(stats_goles$partidos_mas_5_goles)),
  tibble(estadistica = "Marcador más frecuente",          valor = top_marcadores$marcador[1]),
  tibble(estadistica = "Overround promedio (%)",          valor = as.character(round(mean(datos$overround) * 100, 3)))
)

write.csv(tabla_descriptiva,  "tabla_descriptiva_temporada.csv",  row.names = FALSE)
write.csv(goles_por_mes,      "tabla_goles_por_mes.csv",          row.names = FALSE)
write.csv(top_marcadores,     "tabla_marcadores_frecuentes.csv",  row.names = FALSE)
cat("\nArchivos descriptivos guardados.\n")


# -----------------------------------------------------------------------------
# 11. RESUMEN FINAL EN CONSOLA
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("PIPELINE COMPLETO — RESUMEN DE OUTPUTS GENERADOS\n")
cat(strrep("=", 60), "\n\n")

cat("ARCHIVOS .CSV — comparación modelo vs. mercado (validación 24/25):\n")
cat("  tabla_maestra_resultados.csv\n")
cat("  tabla_discrepancia_top10_mayor.csv\n")
cat("  tabla_discrepancia_top10_menor.csv\n")
cat("  tabla_parametros_modelo.csv\n")
cat("  tabla_modelo_top10_aciertos.csv\n")
cat("  tabla_modelo_top10_fallos.csv\n")
cat("  tabla_descriptiva_temporada.csv\n")
cat("  tabla_goles_por_mes.csv\n")
cat("  tabla_marcadores_frecuentes.csv\n\n")

cat("ARCHIVOS .CSV — análisis de mercado (4 temporadas):\n")
cat("  tabla_overround_global.csv\n")
cat("  tabla_overround_temporada.csv\n")
cat("  tabla_overround_pareado.csv\n")
cat("  tabla_overround_robustez_temporal.csv\n")
cat("  tabla_flb_regresiones.csv\n")
cat("  tabla_flb_agregado.csv\n")
cat("  tabla_flb_potencia.csv\n")
cat("  tabla_flb_tendencia.csv\n")
cat("  tabla_retornos_por_banda.csv\n")
cat("  tabla_retornos_extremos.csv\n")
cat("  tabla_correlacion_descriptiva.csv\n\n")

cat("ARCHIVOS .CSV — contrastes principales (Sección 2.1):\n")
cat("  tabla_contraste_primario.csv\n")
cat("  tabla_contrastes_calibracion.csv\n")
cat("  tabla_contraste_sin_ascendidos.csv\n")
cat("  tabla_piso_ruido_rel.csv\n\n")

cat("QUÉ ENCABEZA CADA CAPÍTULO:\n")
cat("  Eficiencia del mercado (FLB) -> tabla_flb_agregado.csv, NO las seis\n")
cat("     regresiones desagregadas.\n")
cat("  Modelo vs. mercado           -> tabla_contraste_primario.csv, no las\n")
cat("     magnitudes sueltas de la tabla maestra.\n")
cat("  Overround pre/cierre         -> leer el contraste global JUNTO con\n")
cat("     tabla_overround_robustez_temporal.csv.\n\n")

cat("GRÁFICOS .PNG (300 dpi):\n")
cat("  grafico_07a_overround_distribucion.png\n")
cat("  grafico_07b_fiabilidad_deciles.png\n")
cat("  grafico_07c_flb_pendientes.png\n")
cat("  grafico_07d_overround_vs_favorito.png\n")
cat("  grafico_07e_overround_temporada.png\n")
cat("  grafico_07f_retornos_por_banda.png\n")
cat("  grafico_07g_potencia_flb.png\n")
cat("  grafico_08a_brier_deciles.png\n")
cat("  grafico_08b_logloss_deciles.png\n")
cat("  grafico_08c_descomposicion_murphy.png\n")
cat("  grafico_08d_brier_por_resultado.png\n")
cat("  grafico_09a_sintesis_modelo_vs_mercado.png\n")
cat("  grafico_09b_parametros_modelo.png\n\n")

cat("Pipeline finalizado.\n")
