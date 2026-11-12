# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 03: Robustez, test de significatividad e intervalo de confianza de ρ
# =============================================================================
# Depende de: Script 02 (resultados_dc.rds)
# Concluye si ρ es significativo y qué modelo usar en el Script 04.
# =============================================================================
#
# DECISIONES DE DISEÑO
#
#   (a) MODELOS ANIDADOS POR CONSTRUCCIÓN.
#       El modelo restringido (Maher) se estima como el modelo DC con ρ FIJADO
#       en 0, reusando la función de verosimilitud del Script 02 mediante un
#       envoltorio, en lugar de escribir una función aparte. Así los dos
#       modelos quedan anidados por construcción y el LRT compara exactamente
#       lo que debe.
#
#       El vector resultado_maher$par tiene el layout de 2n+1 = 51 elementos
#       (m, γ, n−1 ataques, n defensas) que leen los Scripts 04 y 06.
#
#   (b) SE NORMALIZAN LOS PARÁMETROS DEL MODELO FINAL (Sección 4).
#       Misma convención que el Script 02: media(ataques) = 0 y
#       media(defensas) = 0. Se guardan normalizados para que el Script 04 no
#       tenga que repetir la operación ni reporte la media global sesgada por
#       la restricción α₁ = 0.
#
#   (c) LECTURA DE LA PRUEBA DE ROBUSTEZ (Sección 2).
#       Los seis arranques NO convergen al mismo valor: convergen a seis
#       valores distintos. El argumento de robustez no es la coincidencia sino
#       la escala: la dispersión del optimizador es un porcentaje ínfimo del
#       error estándar del parámetro.
#       Por la misma razón no se reporta el promedio de los arranques como si
#       fuera un estimador: son seis aproximaciones numéricas al mismo MLE, no
#       seis observaciones.
#
#   (d) PERFIL DE VEROSIMILITUD DE ρ (Sección 5).
#       El LRT solo dice "no se rechaza ρ = 0", que es un no-hallazgo. El
#       perfil da un INTERVALO DE CONFIANZA para ρ, que es un resultado
#       positivo: acota qué valores de ρ son compatibles con estos datos.
#
#   NOTA SOBRE LOS GRADOS DE LIBERTAD
#       El Script 02 mostró que la parametrización tiene una dirección plana.
#       Eso NO afecta al LRT: ambos modelos comparten exactamente la misma
#       no-identificación, así que la diferencia de parámetros identificados
#       sigue siendo 1 (solo ρ). Además ρ = 0 está en el INTERIOR del espacio
#       de parámetros, no en el borde, de modo que aplica la χ²(1) estándar
#       sin corrección de mezcla.
# =============================================================================
#
# TIEMPO DE EJECUCIÓN
#   Sección 2 (robustez):  6 optimizaciones
#   Sección 3 (LRT):       1 optimización
#   Sección 5 (perfil):    length(grilla_rho) − 1 optimizaciones
#   Con la grilla por defecto son 12 optimizaciones en total (~20-25 minutos).
#   Para una corrida rápida, achicar grilla_rho en la Sección 5.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. CARGAR RESULTADOS DEL SCRIPT 02
# -----------------------------------------------------------------------------

dc <- readRDS("resultados_dc.rds")

resultado_dc         <- dc$resultado_dc
log_verosimilitud_dc <- dc$log_verosimilitud_dc
tau                  <- dc$tau   # la verosimilitud la invoca; sin esta línea el
                                  # script falla al ejecutarse en una sesión nueva de R
datos_train          <- dc$datos_train
equipos              <- dc$equipos
n_equipos            <- dc$n_equipos

loglik_dc <- -resultado_dc$value
rho_hat   <- resultado_dc$par[3]

cat("Log-verosimilitud del modelo con τ:", round(loglik_dc, 2), "\n")
cat("Rho estimado (Script 02):          ", round(rho_hat, 6), "\n")


# -----------------------------------------------------------------------------
# 2. PRUEBA DE ROBUSTEZ DE RHO
# -----------------------------------------------------------------------------
# Se re-optimiza desde distintos valores iniciales de ρ para descartar que la
# estimación sea un artefacto de un óptimo local.
#
# LECTURA CORRECTA: las corridas NO devuelven el mismo valor, sino valores que
# difieren en la cuarta o quinta decimal. Eso es esperable y no es un problema:
# la verosimilitud es casi plana en ρ. Lo que hay que mostrar no es que los
# valores coincidan, sino que su dispersión es despreciable frente a la
# incertidumbre estadística del parámetro. Esa comparación se hace en la
# Sección 3, una vez disponible el error estándar.

cat("\n=== PRUEBA DE ROBUSTEZ DE RHO ===\n")
cat("Se re-optimiza desde distintos valores iniciales de rho.\n\n")

valores_rho_init <- c(-0.5, -0.2, -0.1, 0, 0.1, 0.2)

resultados_robustez <- tibble(
  rho_inicial  = numeric(),
  rho_estimado = numeric(),
  log_verosim  = numeric(),
  convergio    = logical()
)

for (rho_init in valores_rho_init) {

  cat("Rho inicial:", rho_init, "→ optimizando...\n")

  params_test <- c(
    log(mean(datos_train$fthg)),  # media global
    0.1,                          # ventaja local
    rho_init,                     # <-- valor inicial que variamos
    rep(0, n_equipos - 1),        # ataques
    rep(0, n_equipos)             # defensas
  )

  res <- optim(
    par     = params_test,
    fn      = log_verosimilitud_dc,
    datos   = datos_train,
    equipos = equipos,
    method  = "BFGS",
    control = list(maxit = 5000, trace = 0)
  )

  resultados_robustez <- resultados_robustez |>
    add_row(
      rho_inicial  = rho_init,
      rho_estimado = res$par[3],
      log_verosim  = -res$value,
      convergio    = res$convergence == 0
    )
}

cat("\n--- Tabla de resultados ---\n")
resultados_robustez |>
  mutate(rho_estimado = round(rho_estimado, 6),
         log_verosim  = round(log_verosim, 4)) |>
  print()

dispersion_rho <- diff(range(resultados_robustez$rho_estimado))

cat("\n--- Resumen ---\n")
cat("Todos convergieron:                  ", all(resultados_robustez$convergio), "\n")
cat("Dispersión de rho (máx − mín):       ", round(dispersion_rho, 6), "\n")
cat("Dispersión de la log-verosimilitud:  ",
    round(diff(range(resultados_robustez$log_verosim)), 6), "\n")
cat("\nNOTA: el promedio de estas corridas NO es un estimador. Son seis\n")
cat("aproximaciones numéricas al mismo MLE, no seis observaciones. El\n")
cat("estimador puntual es el del Script 02:", round(rho_hat, 6), "\n")
cat("La dispersión se interpreta en la Sección 3, contra el error estándar.\n")

# Gráfico de robustez
rango_y <- range(resultados_robustez$rho_estimado)
resultados_robustez |>
  ggplot(aes(x = rho_inicial, y = rho_estimado)) +
  geom_hline(yintercept = rho_hat, linetype = "dashed", color = "#e74c3c") +
  geom_point(size = 4, color = "#2c3e50") +
  annotate("text",
           x     = min(valores_rho_init),
           y     = rho_hat + 0.15 * diff(rango_y),
           label = paste("Estimación del Script 02:", round(rho_hat, 5)),
           color = "#e74c3c", hjust = 0, size = 3.5) +
  labs(
    title    = "Robustez del parámetro ρ frente al valor inicial",
    subtitle = "Poisson con corrección τ | Premier League 2021/22 - 2023/24",
    x        = "Valor inicial de rho",
    y        = "Rho estimado"
  ) +
  theme_minimal()


# -----------------------------------------------------------------------------
# 3. TEST DE RAZÓN DE VEROSIMILITUDES (LRT)
# -----------------------------------------------------------------------------
# H0: ρ = 0  (modelo de Maher, 1982)
# H1: ρ ≠ 0  (modelo con corrección τ de Dixon & Coles, 1997)
# Estadístico: LR = 2 (logLik_τ − logLik_Maher) ~ χ²(1) bajo H0
#
# El modelo restringido se estima FIJANDO ρ en 0 dentro de la misma función de
# verosimilitud del Script 02. No hay una segunda función escrita a mano.

cat("\n=== TEST DE RAZÓN DE VEROSIMILITUDES ===\n")

# Envoltorio: fija ρ en un valor dado y deja que optim() mueva el resto.
# params_sin_rho = (media_global, ventaja_local, n−1 ataques, n defensas)
log_verosim_rho_fijo <- function(params_sin_rho, datos, equipos, rho_fijo) {
  log_verosimilitud_dc(
    c(params_sin_rho[1:2], rho_fijo, params_sin_rho[-(1:2)]),
    datos, equipos
  )
}

estimar_con_rho_fijo <- function(rho_fijo) {
  optim(
    par      = c(log(mean(datos_train$fthg)), 0.1,
                 rep(0, n_equipos - 1), rep(0, n_equipos)),
    fn       = log_verosim_rho_fijo,
    datos    = datos_train,
    equipos  = equipos,
    rho_fijo = rho_fijo,
    method   = "BFGS",
    control  = list(maxit = 5000, trace = 0)
  )
}

cat("Estimando el modelo restringido (Maher, rho = 0)...\n")

resultado_maher <- estimar_con_rho_fijo(0)
loglik_maher    <- -resultado_maher$value

cat("¿Convergió?", if (resultado_maher$convergence == 0) "SÍ ✓" else "NO - revisar", "\n")
cat("Log-verosimilitud Maher (rho = 0):", round(loglik_maher, 4), "\n")
cat("Log-verosimilitud con τ:          ", round(loglik_dc,    4), "\n")

# --- Chequeo de verificación 1: el modelo restringido no puede ajustar mejor --
# Por anidamiento, logLik_τ ≥ logLik_Maher siempre. Una violación indicaría que
# la optimización del modelo completo se detuvo antes de tiempo.
if (loglik_dc < loglik_maher) {
  warning("El modelo restringido ajusta MEJOR que el completo (LR < 0). ",
          "La optimización del Script 02 no convergió al óptimo. ",
          "Revisar antes de interpretar el test.")
} else {
  cat("Chequeo de anidamiento: logLik_τ >= logLik_Maher ✓\n")
}

# --- Estadístico LR y p-valor ---
LR      <- 2 * (loglik_dc - loglik_maher)
p_valor <- pchisq(LR, df = 1, lower.tail = FALSE)

cat("\n--- Resultado del test ---\n")
cat("Estadístico LR:     ", round(LR,      4), "\n")
cat("Grados de libertad:  1\n")
cat("P-valor:            ", round(p_valor, 4), "\n")

cat("\n--- Interpretación ---\n")
if (p_valor < 0.05) {
  cat("Rho es estadísticamente significativo (p =", round(p_valor, 4), ")\n")
  cat("  → Se justifica la corrección τ.\n")
} else {
  cat("Rho NO es estadísticamente significativo (p =", round(p_valor, 4), ")\n")
  cat("  → El modelo de Maher (sin corrección τ) es suficiente para estos datos.\n")
}

# --- Error estándar aproximado de rho y lectura de la robustez ----------------
# Bajo aproximación cuadrática del perfil de verosimilitud,
#   LR ≈ (rho_hat / EE)²  →  EE ≈ |rho_hat| / sqrt(LR)
# Este EE usa las dos estimaciones que ya tenemos (rho libre y rho = 0), sin
# costo adicional. La Sección 5 lo valida contra el perfil completo.

ee_rho_aprox <- abs(rho_hat) / sqrt(LR)

cat("\n--- Escala de la dispersión del optimizador ---\n")
cat("Error estándar aproximado de rho:", round(ee_rho_aprox, 4), "\n")
cat(sprintf("Dispersión entre corridas (Sección 2): %.6f = %.2f%% de un EE\n",
            dispersion_rho, dispersion_rho / ee_rho_aprox * 100))
cat("-> La dispersión numérica del optimizador es irrelevante frente a la\n")
cat("   incertidumbre estadística del parámetro. ESE es el argumento de\n")
cat("   robustez, no que las corridas devuelvan el mismo número.\n")


# -----------------------------------------------------------------------------
# 4. PARÁMETROS DEL MODELO FINAL, NORMALIZADOS
# -----------------------------------------------------------------------------
# Misma convención de identificación que el Script 02: media(ataques) = 0 y
# media(defensas) = 0. Se normaliza acá para que el Script 04 reciba los
# parámetros ya en escala interpretable y no reporte la media global sesgada
# por la restricción α₁ = 0 (que la vuelve la tasa del primer equipo
# alfabético en lugar de la de un equipo promedio).

p_m             <- resultado_maher$par
media_global_m  <- p_m[1]
ventaja_local_m <- p_m[2]
ataques_m       <- c(0, p_m[3:(n_equipos + 1)])
defensas_m      <- p_m[(n_equipos + 2):(2 * n_equipos + 1)]
names(ataques_m)  <- equipos
names(defensas_m) <- equipos

c_ataque_m  <- mean(ataques_m)
c_defensa_m <- mean(defensas_m)

ataques_m      <- ataques_m  - c_ataque_m
defensas_m     <- defensas_m - c_defensa_m
media_global_m <- media_global_m + c_ataque_m + c_defensa_m

cat("\n--- Parámetros del modelo final (Maher), normalizados ---\n")
cat("Media global (log):", round(media_global_m, 4), "\n")
cat("  → equipo promedio de local:    ",
    round(exp(media_global_m + ventaja_local_m), 3), "goles\n")
cat("  → equipo promedio de visitante:",
    round(exp(media_global_m), 3), "goles\n")
cat("Ventaja local (log):", round(ventaja_local_m, 4),
    "→ factor multiplicativo:", round(exp(ventaja_local_m), 3), "\n")
cat("Medias tras normalizar → ataques:", round(mean(ataques_m), 10),
    "| defensas:", round(mean(defensas_m), 10), "(ambas deben ser 0)\n")

# --- Chequeo de verificación 2: medias ajustadas vs observadas ---------------
# Con link logarítmico e intercepto, las ecuaciones de score fuerzan la
# coincidencia. Una discrepancia apreciable indica error de indexación.
lambdas_aj <- exp(media_global_m + ataques_m[datos_train$home_team] +
                  defensas_m[datos_train$away_team] + ventaja_local_m)
mus_aj     <- exp(media_global_m + ataques_m[datos_train$away_team] +
                  defensas_m[datos_train$home_team])

cat(sprintf("\nMedias ajustadas vs observadas → local: %.4f / %.4f | visitante: %.4f / %.4f\n",
            mean(lambdas_aj), mean(datos_train$fthg),
            mean(mus_aj),     mean(datos_train$ftag)))

if (max(abs(c(mean(lambdas_aj) - mean(datos_train$fthg),
              mean(mus_aj)     - mean(datos_train$ftag)))) > 0.02) {
  warning("Las medias ajustadas no reproducen las observadas. ",
          "Revisar los índices de extracción de ataques y defensas.")
} else {
  cat("Chequeo de medias ajustadas ✓\n")
}


# -----------------------------------------------------------------------------
# 5. PERFIL DE VEROSIMILITUD DE ρ E INTERVALO DE CONFIANZA
# -----------------------------------------------------------------------------
# El LRT de la Sección 3 solo concluye "no se rechaza ρ = 0", que es un
# no-hallazgo. El perfil de verosimilitud da un resultado positivo: acota qué
# valores de ρ son compatibles con estos datos.
#
# Para cada ρ de la grilla se fija ese valor y se re-optimizan los 51
# parámetros restantes. La curva resultante es la verosimilitud perfilada.
# El intervalo de confianza al 95% son los ρ con 2(logLik_max − logLik_perfil)
# menor que el cuantil 0,95 de una χ²(1) = 3,841.
#
# COSTO: una optimización por punto de la grilla. Achicar si hace falta.

grilla_rho <- c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15)

cat("\n=== PERFIL DE VEROSIMILITUD DE RHO ===\n")
cat("Optimizando", length(grilla_rho), "puntos de la grilla...\n")

perfil_rho <- map_dfr(grilla_rho, function(r) {
  cat("  rho fijo =", r, "\n")
  # Reutiliza el ajuste de la Sección 3 en lugar de repetirlo
  res <- if (r == 0) resultado_maher else estimar_con_rho_fijo(r)
  tibble(rho = r, loglik_perfil = -res$value, convergio = res$convergence == 0)
})

perfil_rho <- perfil_rho |>
  mutate(LR_perfil = 2 * (loglik_dc - loglik_perfil))

cat("\n--- Perfil de verosimilitud ---\n")
print(perfil_rho |> mutate(across(where(is.numeric), \(x) round(x, 4))))

if (!all(perfil_rho$convergio)) {
  warning("Algún punto de la grilla no convergió. Revisar antes de interpretar el IC.")
}

# --- Intervalo de confianza al 95% por interpolación sobre el perfil ---------
umbral_chi2 <- qchisq(0.95, df = 1)   # 3.841

rama_izq <- perfil_rho |> filter(rho <= rho_hat)
rama_der <- perfil_rho |> filter(rho >= rho_hat)

ic_inf <- if (max(rama_izq$LR_perfil) >= umbral_chi2) {
  approx(x = rama_izq$LR_perfil, y = rama_izq$rho, xout = umbral_chi2)$y
} else NA_real_

ic_sup <- if (max(rama_der$LR_perfil) >= umbral_chi2) {
  approx(x = rama_der$LR_perfil, y = rama_der$rho, xout = umbral_chi2)$y
} else NA_real_

cat("\n--- Intervalo de confianza al 95% para rho ---\n")
if (is.na(ic_inf) || is.na(ic_sup)) {
  cat("La grilla no cubre el umbral chi2 en alguna rama. Ampliar grilla_rho:\n")
  cat("  LR máximo alcanzado - rama izquierda:", round(max(rama_izq$LR_perfil), 3),
      "| rama derecha:", round(max(rama_der$LR_perfil), 3),
      "| umbral:", round(umbral_chi2, 3), "\n")
} else {
  cat(sprintf("  rho_hat = %.5f   IC95 = (%.4f , %.4f)   amplitud = %.4f\n",
              rho_hat, ic_inf, ic_sup, ic_sup - ic_inf))
}

# --- Validación de la aproximación cuadrática usada en la Sección 3 ----------
# Si el perfil es aproximadamente cuadrático, el IC del perfil y el IC de Wald
# construido con ee_rho_aprox deben coincidir de cerca.
cat("\n--- Contraste con la aproximación cuadrática (Sección 3) ---\n")
cat(sprintf("  IC95 de Wald aproximado: (%.4f , %.4f)\n",
            rho_hat - qnorm(0.975) * ee_rho_aprox,
            rho_hat + qnorm(0.975) * ee_rho_aprox))
cat("  Si ambos intervalos son parecidos, el perfil es casi cuadrático y el EE\n")
cat("  aproximado de la Sección 3 es utilizable para reportar.\n")

# --- Gráfico del perfil ------------------------------------------------------
perfil_rho |>
  ggplot(aes(x = rho, y = LR_perfil)) +
  geom_line(color = "#2c3e50") +
  geom_point(size = 2.5, color = "#2c3e50") +
  geom_hline(yintercept = umbral_chi2, linetype = "dashed", color = "#e74c3c") +
  geom_vline(xintercept = rho_hat, linetype = "dotted", color = "#2c3e50") +
  annotate("text", x = min(grilla_rho), y = umbral_chi2,
           label = "χ²(1) al 95% = 3,841", color = "#e74c3c",
           hjust = 0, vjust = -0.5, size = 3.5) +
  labs(
    title    = "Verosimilitud perfilada del parámetro ρ",
    subtitle = "Los ρ por debajo de la línea son compatibles con los datos al 95%",
    x        = "ρ",
    y        = "Estadístico LR = 2 (logLik máx − logLik perfilada)"
  ) +
  theme_minimal()


# -----------------------------------------------------------------------------
# 6. GUARDAR PARA EL SCRIPT 04
# -----------------------------------------------------------------------------
# resultado_maher$par tiene el layout de 2n+1 elementos que leen los
# Scripts 04 y 06. Los parámetros normalizados van aparte.

# log_verosimilitud_maher se guarda porque el Script 04 lo carga (aunque no lo
# invoque). Es un envoltorio sobre la función del Script 02, no una
# reimplementación.
log_verosimilitud_maher <- function(params, datos, equipos) {
  log_verosim_rho_fijo(params, datos, equipos, rho_fijo = 0)
}

saveRDS(
  list(
    log_verosimilitud_maher = log_verosimilitud_maher,
    log_verosim_rho_fijo    = log_verosim_rho_fijo,
    loglik_maher            = loglik_maher,
    loglik_dc               = loglik_dc,
    LR                      = LR,
    p_valor                 = p_valor,
    resultado_maher         = resultado_maher,
    # --- inferencia sobre rho ---
    rho_hat                 = rho_hat,
    ee_rho_aprox            = ee_rho_aprox,
    perfil_rho              = perfil_rho,
    ic95_rho                = c(inf = ic_inf, sup = ic_sup),
    # --- parámetros normalizados: media(ataques) = 0 y media(defensas) = 0 ---
    media_global_norm       = media_global_m,
    ventaja_local           = ventaja_local_m,
    ataques_norm            = ataques_m,
    defensas_norm           = defensas_m,
    c_ataque                = c_ataque_m,
    c_defensa               = c_defensa_m,
    # --- robustez ---
    resultados_robustez     = resultados_robustez
  ),
  "resultados_lrt.rds"
)

cat("\nArchivo 'resultados_lrt.rds' guardado → usar en Script 04.\n")
