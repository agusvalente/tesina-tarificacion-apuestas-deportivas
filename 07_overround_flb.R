# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 07: Overround y Favourite-Longshot Bias
# =============================================================================
# Depende de: Script 01 (datos_mercado.rds)
#             Script 06 (resultados_modelo_final.rds)
#
# Objetivos:
#   1. Calcular el overround (loading actuarial) en las DOS series de cuotas
#      Bet365 —pre-cierre y cierre— sobre las cuatro temporadas
#   2. Comparar ambas series de forma pareada (mismos partidos)
#   3. Normalizar las probabilidades implícitas de mercado
#   4. Contrastar formalmente el favourite-longshot bias mediante regresión
#      logística de calibración, por categoría y por serie, y mediante un test
#      agregado con errores estándar agrupados por partido
#   5. Cuantificar la potencia del test y la precisión de las estimaciones
#   6. Conservar el diagrama de fiabilidad por deciles como visualización
#   7. Replicar el test de retornos por banda de precio de Cain, Law y Peel (2000)
#   8. Dejar preparado el objeto de validación 2024/25 para los Scripts 08 y 09
#
# -----------------------------------------------------------------------------
# DECISIONES DE DISEÑO
#
#   (a) DOS SERIES DE CUOTAS. Las cuotas b365h/d/a NO son de cierre sino las
#       recogidas el viernes previo. Se reportan las dos series: pre-cierre
#       (serie primaria, es el precio real que enfrenta el apostador minorista
#       durante varios días, coherente con la lectura del overround como
#       loading actuarial) y cierre (robustez).
#       AMBAS SE REPORTAN SIEMPRE, den lo que den. Elegir después la que arroje
#       el resultado más conveniente sería data dredging, el mismo problema que
#       se evita al pre-especificar LRT/AIC/BIC en la selección del modelo.
#
#   (b) VENTANA TEMPORAL. El análisis de mercado usa las CUATRO temporadas
#       (≈1.520 partidos): no involucra al modelo, solo cuotas y resultados
#       observados, por lo que no necesita ser fuera de muestra. La restricción
#       de usar únicamente 2024/25 aplica a la comparación modelo vs. mercado
#       (Script 08), no al estudio del mercado en sí.
#
#   (c) TEST DE FLB. La inferencia la hace una regresión logística de
#       calibración sobre los datos individuales. Agrupar las observaciones en
#       10 deciles y correlacionar los 10 puntos descartaría casi toda la
#       información, no daría error estándar ni test formal, y las filas no son
#       independientes (cada partido aporta 3 categorías que suman 1). Los
#       deciles quedan como gráfico, y la correlación como estadístico
#       descriptivo.
#
#   (d) DIAGNÓSTICO DE VALORES EXTREMOS DEL OVERROUND (Sección 3.1).
#       La serie pre-cierre registra un máximo de 16,7% contra una media de
#       5,43% y un desvío de 0,65: más de 17 desvíos, mientras que la serie de
#       cierre no supera el 8,06%. Un margen así no es creíble en Premier
#       League y apunta a un error de carga del proveedor. Se identifica, se
#       reporta, y el test pareado se corre con y sin esos partidos
#       (Sección 4.2). No se elimina nada del análisis principal: se muestra
#       si la conclusión depende o no del dato dudoso.
#
#   (e) TEST AGREGADO CON ERRORES ESTÁNDAR AGRUPADOS (Sección 7.2).
#       Las seis regresiones por categoría y serie dan seis respuestas
#       parciales y correlacionadas entre sí, lo que impide sintetizarlas
#       mirando p-valores uno por uno. Por eso se estima además una regresión
#       única con interceptos por categoría y pendiente común, con varianza
#       agrupada a nivel partido. Devuelve un único contraste con inferencia
#       correcta.
#
#   (f) LECTURA POR PRECISIÓN DEL INTERVALO (Sección 7.4). Se reporta
#       explícitamente qué valores de b descarta y cuáles no descarta cada
#       intervalo de confianza. Es preferible a evaluar la potencia en el valor
#       estimado —práctica desaconsejada, porque la potencia observada es
#       función monótona del p-valor y no aporta información—. La tabla de
#       potencia se mantiene porque evalúa valores PRE-ESPECIFICADOS.
#
#   (g) REFERENCIA DE RETORNO JUSTO (Sección 9). La línea de referencia del
#       gráfico de Cain no está en cero: con overround ninguna apuesta puede
#       rendir cero. Bajo calibración perfecta y margen uniforme toda apuesta
#       unitaria rinde 1/(1+OR) − 1, independientemente de la cuota. Esa es la
#       referencia correcta y hace visible el sesgo de inmediato.
#
#   (h) CRITERIOS DE LECTURA. La lectura actuarial de la Sección 4 considera
#       tanto el test t como el de Wilcoxon, y el chequeo de coherencia de la
#       Sección 9 compara dirección y magnitud, no banderas de significación,
#       para que una diferencia marginal no se lea como una contradicción.
#
#   (i) ERRORES ESTÁNDAR AGRUPADOS EN TODA LA SECCIÓN 9. Como en la Sección
#       7.2, la tabla de retornos por banda, el contraste entre bandas extremas
#       y el test de tendencia agrupan por partido.
#
#       El caso crítico es el contraste entre extremos: las dos bandas están
#       pobladas por los MISMOS partidos —un partido con un favorito a 1,40
#       (i ≈ 0,71) tiene por construcción un visitante a 7,00 (i ≈ 0,14)— y
#       dentro del partido los retornos están fuertemente correlacionados en
#       sentido negativo. Como Var(R_alta − R_baja) = Var_a + Var_b − 2·Cov
#       con Cov < 0, el error estándar ingenuo SUBESTIMA la incertidumbre y el
#       contraste sería anti-conservador. Por eso se estima sobre las
#       observaciones individuales con varianza agrupada por partido.
#
#   (j) LECTURA EN ZONA DE BORDE (Sección 7.4.2). Dicotomizar en 0,05 llevaría,
#       con p = 0,0451 en una regresión individual y p = 0,0513 en el test
#       agregado, a afirmar que el resultado "no sobrevive a la agregación",
#       lo que es falso: son en la práctica el mismo número. Hay una rama
#       explícita para el caso en que el agregado queda en la zona 0,05–0,10
#       con la pendiente del lado del FLB, coherente con el criterio de la
#       Sección 9: comparar dirección y magnitud, no banderas de
#       significación.
#
#   (k) GRILLA DE POTENCIA. Incluye 1,15 y 1,25 además de 1,10 / 1,20 / 1,30 /
#       1,50. Con solo los valores gruesos el efecto mínimo detectable de H y
#       A quedaría en 1,30, cuando la potencia en 1,20 ya es ≈ 0,79: reportar
#       "solo detecta b ≥ 1,30" subestimaría el diseño.
#
#   (l) NOTAS DE LIMITACIÓN. (1) Los intervalos de Wilson de la Sección 8 no
#       agrupan por partido, por lo que el recuento de deciles fuera del IC es
#       levemente anti-conservador: se declara y el bloque queda como
#       estrictamente descriptivo. (2) La correlación π–sesgo por decil CAMBIA
#       DE SIGNO entre las cuatro temporadas y la temporada de validación
#       2024/25, lo que condiciona la lectura de los Scripts 08 y 09: se
#       reporta explícitamente en la Sección 10.
#
#   (m) COHERENCIA DE LA LECTURA ACTUARIAL (Sección 4.5). La lectura declara
#       primero el ALCANCE TEMPORAL del resultado —si sobrevive o no a excluir
#       2021/22, que establece la Sección 4.4— y recién después la dirección,
#       para que ninguna rama decida sin mirar la evidencia que el propio
#       script ya produjo.
#
# -----------------------------------------------------------------------------
# NOTA SOBRE EL MÉTODO DE NORMALIZACIÓN (para el capítulo de limitaciones)
#
#   Se usa normalización proporcional (dividir cada probabilidad implícita
#   cruda por la suma de las tres). Shin (1993) propone una alternativa que
#   descuenta el margen atribuible a la presencia de agentes informados y que,
#   en principio, podría inducir o enmascarar parte del sesgo medido.
#
#   No se implementa acá, y se declara como limitación. El respaldo empírico
#   de que la decisión no es determinante lo aporta el propio antecedente que
#   se cita en el marco teórico: Cain, Law y Peel (2000, nota al pie 4)
#   reportan que las probabilidades ajustadas por el modelo de Shin y el
#   precio crudo tienen R² ≈ 0,9994 en una regresión lineal
#   (π = 0,027 + 1,026 · SHINP), es decir, son prácticamente colineales.
#   Bajo esa evidencia, el método de normalización es poco probable que altere
#   las conclusiones de forma sustantiva.
#
# -----------------------------------------------------------------------------
# ORDEN DE EJECUCIÓN
#   El script corre de arriba a abajo sin saltos. Cada sección usa solo
#   objetos definidos antes. Tiempo aproximado: 1 a 2 minutos, casi todo en
#   la simulación de potencia de la Sección 7.4.
# =============================================================================

library(tidyverse)


# -----------------------------------------------------------------------------
# 1. CARGAR DATOS
# -----------------------------------------------------------------------------
# Dos insumos con dos ventanas temporales distintas (ver nota (b)):
#   datos_mercado  -> 4 temporadas, análisis de mercado
#   modelo         -> 380 partidos de 24/25 con probabilidades del modelo

datos_mercado <- readRDS("datos_mercado.rds")
modelo        <- readRDS("resultados_modelo_final.rds")

datos_validacion <- modelo$resultados_fin

# Salvaguarda: los Scripts 02 a 06 pueden haber descartado columnas por el
# camino. Si las cuotas de cierre no llegaron hasta acá, se recuperan desde
# datos_mercado uniendo por partido. No se toca nada más del objeto.
if (!all(c("b365ch", "b365cd", "b365ca") %in% names(datos_validacion))) {
  cat("Las cuotas de cierre no están en resultados_fin: se recuperan desde datos_mercado.\n")
  datos_validacion <- datos_validacion |>
    left_join(
      datos_mercado |>
        select(date, home_team, away_team, b365ch, b365cd, b365ca),
      by = c("date", "home_team", "away_team")
    )
  stopifnot(nrow(datos_validacion) == nrow(modelo$resultados_fin))
}

cat("=== INSUMOS ===\n")
cat("Análisis de mercado (4 temporadas):", nrow(datos_mercado), "partidos\n")
cat("Validación 24/25 (modelo):         ", nrow(datos_validacion), "partidos\n")

# Cobertura de las cuotas de cierre
n_sin_cierre <- sum(is.na(datos_mercado$b365ch) |
                    is.na(datos_mercado$b365cd) |
                    is.na(datos_mercado$b365ca))
cat("Partidos sin cuota de cierre:      ", n_sin_cierre,
    sprintf("(%.2f%%)\n", n_sin_cierre / nrow(datos_mercado) * 100))


# -----------------------------------------------------------------------------
# 2. OVERROUND Y PROBABILIDADES IMPLÍCITAS NORMALIZADAS
# -----------------------------------------------------------------------------
# Las probabilidades implícitas crudas son simplemente 1/cuota.
# Si el mercado fuera justo, sumarían exactamente 1.
# El exceso sobre 1 es el overround: el margen de la casa de apuestas,
# equivalente al loading actuarial en tarificación de seguros.
#
# La normalización proporcional divide cada probabilidad cruda por la suma de
# las tres. Se aplica ANTES del análisis de FLB: sin ella, el margen inflaría
# todas las probabilidades implícitas y generaría un sesgo sistemático espurio
# (las tres categorías aparecerían "sobrestimadas" solo por el recargo).
#
# Convención de nombres:
#   pi_H,  pi_D,  pi_A   -> serie PRE-CIERRE (primaria)
#   piC_H, piC_D, piC_A  -> serie de CIERRE  (robustez)

agregar_probabilidades <- function(df) {
  df |>
    mutate(
      # --- Serie pre-cierre ---
      pi_H_cruda = 1 / b365h,
      pi_D_cruda = 1 / b365d,
      pi_A_cruda = 1 / b365a,
      suma_pi    = pi_H_cruda + pi_D_cruda + pi_A_cruda,
      overround  = suma_pi - 1,
      pi_H       = pi_H_cruda / suma_pi,
      pi_D       = pi_D_cruda / suma_pi,
      pi_A       = pi_A_cruda / suma_pi,

      # --- Serie de cierre ---
      piC_H_cruda = 1 / b365ch,
      piC_D_cruda = 1 / b365cd,
      piC_A_cruda = 1 / b365ca,
      suma_piC    = piC_H_cruda + piC_D_cruda + piC_A_cruda,
      overround_c = suma_piC - 1,
      piC_H       = piC_H_cruda / suma_piC,
      piC_D       = piC_D_cruda / suma_piC,
      piC_A       = piC_A_cruda / suma_piC
    )
}

datos_mercado    <- agregar_probabilidades(datos_mercado)
datos_validacion <- agregar_probabilidades(datos_validacion)

# Verificación: las probabilidades normalizadas deben sumar exactamente 1
cat("\nVerificación — suma de probabilidades normalizadas (debe ser 1):\n")
datos_mercado |>
  summarise(
    pre_cierre_min = min(round(pi_H + pi_D + pi_A, 6)),
    pre_cierre_max = max(round(pi_H + pi_D + pi_A, 6)),
    cierre_min     = min(round(piC_H + piC_D + piC_A, 6), na.rm = TRUE),
    cierre_max     = max(round(piC_H + piC_D + piC_A, 6), na.rm = TRUE)
  ) |>
  print()


# -----------------------------------------------------------------------------
# 3. OVERROUND: DESCRIPTIVOS POR SERIE Y POR TEMPORADA
# -----------------------------------------------------------------------------

resumen_overround <- function(x, etiqueta) {
  x <- x[!is.na(x)]
  tibble(
    serie   = etiqueta,
    n       = length(x),
    media   = mean(x)   * 100,
    mediana = median(x) * 100,
    desvio  = sd(x)     * 100,
    minimo  = min(x)    * 100,
    maximo  = max(x)    * 100
  )
}

overround_global <- bind_rows(
  resumen_overround(datos_mercado$overround,   "Pre-cierre"),
  resumen_overround(datos_mercado$overround_c, "Cierre")
) |>
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n=== OVERROUND GLOBAL (4 temporadas, en %) ===\n")
print(overround_global)

overround_temporada <- datos_mercado |>
  summarise(
    n              = n(),
    pre_cierre     = round(mean(overround,   na.rm = TRUE) * 100, 3),
    cierre         = round(mean(overround_c, na.rm = TRUE) * 100, 3),
    .by = temporada
  ) |>
  mutate(diferencia = round(cierre - pre_cierre, 3)) |>
  arrange(temporada)

cat("\n=== OVERROUND MEDIO POR TEMPORADA (en %) ===\n")
print(overround_temporada)

# -----------------------------------------------------------------------------
# 3.1 DIAGNÓSTICO DE VALORES EXTREMOS DEL OVERROUND
# -----------------------------------------------------------------------------
# Un overround muy por encima del resto de la distribución no es un dato
# interesante: es casi siempre un error de carga del proveedor (una cuota con
# la coma corrida). Conviene identificarlo antes de correr cualquier test,
# porque una sola observación de ese tipo puede dominar la varianza.
#
# No afecta al análisis de FLB: la normalización proporcional divide el margen,
# de modo que las probabilidades normalizadas de ese partido siguen sumando 1.
# Sí afecta a los descriptivos del overround y al test pareado de la Sección 4.

cat("\n", strrep("=", 78), "\n", sep = "")
cat("DIAGNÓSTICO DE VALORES EXTREMOS DEL OVERROUND\n")
cat(strrep("=", 78), "\n", sep = "")

umbral_or <- 0.10   # 10%: por encima de eso el margen es implausible

extremos <- datos_mercado |>
  mutate(
    or_pre = overround   * 100,
    or_cie = overround_c * 100,
    dif_or = (overround_c - overround) * 100
  ) |>
  filter(or_pre > umbral_or * 100 | or_cie > umbral_or * 100 | abs(dif_or) > 3) |>
  select(temporada, date, home_team, away_team,
         b365h, b365d, b365a, or_pre,
         b365ch, b365cd, b365ca, or_cie, dif_or) |>
  arrange(desc(or_pre))

cat(sprintf("\n--- Partidos con overround > %.0f%% en alguna serie, o diferencia > 3 pp ---\n",
            umbral_or * 100))

if (nrow(extremos) == 0) {
  cat("Ninguno. La distribución del overround no presenta valores anómalos.\n")
} else {
  extremos |>
    mutate(across(where(is.numeric), ~ round(., 3))) |>
    print(n = Inf, width = Inf)
  cat("\nRevisar estas cuotas contra el sitio de football-data.co.uk.\n")
  cat("Si alguna es implausible, es un error de carga del proveedor y debe\n")
  cat("declararse en el texto junto con el tratamiento que se le dio.\n")
}

# Comparación de descriptivos con y sin los valores extremos
cat("\n--- Descriptivos del overround con y sin valores extremos ---\n")
bind_rows(
  overround_global |> mutate(version = "Muestra completa"),
  bind_rows(
    resumen_overround(
      datos_mercado$overround[datos_mercado$overround < umbral_or], "Pre-cierre"),
    resumen_overround(
      datos_mercado$overround_c[datos_mercado$overround_c < umbral_or], "Cierre")
  ) |>
    mutate(across(where(is.numeric), ~ round(., 3)), version = "Sin extremos")
) |>
  select(version, serie, n, media, mediana, desvio, minimo, maximo) |>
  arrange(serie, version) |>
  print(width = Inf)

cat("\nSi la MEDIANA apenas cambia y la MEDIA o el DESVÍO sí, se confirma que\n")
cat("el problema es puntual y no una característica de la distribución.\n")



# -----------------------------------------------------------------------------
# 4. COMPARACIÓN PAREADA ENTRE SERIES
# -----------------------------------------------------------------------------
# 4.1 TESTS PAREADOS SOBRE LA MUESTRA COMPLETA
# -----------------------------------------------------------------------------
# No es un control técnico: es un resultado con contenido propio.
# Muestra cómo evoluciona el recargo a medida que se resuelve la incertidumbre
# entre el viernes previo y el momento del cierre, que es exactamente la
# lectura actuarial del overround como loading sobre la prima pura.
#
# Al ser los MISMOS partidos, corresponde un test pareado. Se reportan el
# t pareado (paramétrico) y el de Wilcoxon (no paramétrico, robusto a
# asimetría en la distribución de las diferencias).

pareado <- datos_mercado |>
  filter(!is.na(overround), !is.na(overround_c)) |>
  mutate(dif = overround_c - overround)

cat("\n=== COMPARACIÓN PAREADA: CIERRE − PRE-CIERRE ===\n")
cat("Partidos con ambas series:", nrow(pareado), "\n")
cat(sprintf("Diferencia media:   %+.4f pp\n", mean(pareado$dif) * 100))
cat(sprintf("Diferencia mediana: %+.4f pp\n", median(pareado$dif) * 100))
cat(sprintf("Desvío de la dif.:   %.4f pp\n", sd(pareado$dif) * 100))

t_pareado <- t.test(pareado$overround_c, pareado$overround, paired = TRUE)
w_pareado <- wilcox.test(pareado$overround_c, pareado$overround,
                         paired = TRUE, exact = FALSE)

cat("\n--- Test t pareado ---\n")
cat(sprintf("  t = %.4f | gl = %d | p = %.6f\n",
            t_pareado$statistic, t_pareado$parameter, t_pareado$p.value))
cat(sprintf("  IC 95%% de la diferencia media: [%+.4f ; %+.4f] pp\n",
            t_pareado$conf.int[1] * 100, t_pareado$conf.int[2] * 100))

cat("\n--- Test de Wilcoxon pareado ---\n")
cat(sprintf("  V = %.0f | p = %.6f\n", w_pareado$statistic, w_pareado$p.value))

cat("\nInterpretación: diferencia negativa = el margen se COMPRIME hacia el cierre\n")
cat("                (el recargo cae a medida que se resuelve la incertidumbre).\n")

# Test pareado también por temporada
pareado_temporada <- pareado |>
  summarise(
    n       = n(),
    dif_pp  = round(mean(dif) * 100, 4),
    p_valor = round(t.test(overround_c, overround, paired = TRUE)$p.value, 6),
    .by = temporada
  ) |>
  arrange(temporada)

# -----------------------------------------------------------------------------
# 4.2 DESAGREGACIÓN POR TEMPORADA
# -----------------------------------------------------------------------------
cat("\n--- Comparación pareada por temporada ---\n")
print(pareado_temporada)



# -----------------------------------------------------------------------------
# 4.3 ROBUSTEZ (I): TEST PAREADO SIN LOS VALORES EXTREMOS
# -----------------------------------------------------------------------------
# El test t usa la media y el desvío de las diferencias, de modo que una
# observación extrema lo afecta doblemente. El Wilcoxon trabaja con rangos y
# es inmune. Cuando ambos discrepan, la primera hipótesis a descartar es
# justamente la presencia de un valor anómalo.
#
# No se elimina nada del análisis principal: se reporta el test con y sin,
# que es la forma correcta de tratar un dato sospechoso.

pareado_robusto <- pareado |>
  filter(overround < umbral_or, overround_c < umbral_or)

n_excluidos  <- nrow(pareado) - nrow(pareado_robusto)
hay_extremos <- n_excluidos > 0

cat(sprintf("\n--- Test pareado excluyendo %d partido(s) con overround > %.0f%% ---\n",
            n_excluidos, umbral_or * 100))

if (hay_extremos) {

  t_rob <- t.test(pareado_robusto$overround_c, pareado_robusto$overround,
                  paired = TRUE)
  w_rob <- wilcox.test(pareado_robusto$overround_c, pareado_robusto$overround,
                       paired = TRUE, exact = FALSE)

  comparacion_robustez <- tibble(
    version    = c("Muestra completa", "Sin valores extremos"),
    n          = c(nrow(pareado), nrow(pareado_robusto)),
    dif_pp     = c(mean(pareado$dif), mean(pareado_robusto$dif)) * 100,
    sd_pp      = c(sd(pareado$dif),   sd(pareado_robusto$dif))   * 100,
    t          = c(t_pareado$statistic, t_rob$statistic),
    p_t        = c(t_pareado$p.value,   t_rob$p.value),
    p_wilcoxon = c(w_pareado$p.value,   w_rob$p.value)
  )

  comparacion_robustez |>
    mutate(across(where(is.numeric), ~ round(., 5))) |>
    print(width = Inf)

  cambia_conclusion <- (t_pareado$p.value < 0.05) != (t_rob$p.value < 0.05)

  if (cambia_conclusion) {
    cat("\nLa conclusión del test t CAMBIA al excluir el valor extremo.\n")
    cat("Hay que reportar ambas versiones y explicar el dato sospechoso.\n")
    cat("El Wilcoxon, al trabajar con rangos, ya era inmune al problema:\n")
    cat("por eso discrepaba del t en la muestra completa. La versión sin\n")
    cat("extremos y el Wilcoxon coinciden, lo que resuelve la discrepancia.\n")
  } else {
    cat("\nLa conclusión del test t NO cambia al excluir el valor extremo.\n")
    cat("Si el t y el Wilcoxon siguen discrepando, el origen es la asimetría\n")
    cat("general de la distribución de diferencias, no un dato puntual.\n")
  }

} else {
  cambia_conclusion <- FALSE
  cat("No hay partidos por encima del umbral: nada que excluir.\n")
}



# -----------------------------------------------------------------------------
# 4.4 ROBUSTEZ (II): ¿EL EFECTO DEPENDE DE UNA SOLA TEMPORADA?
# -----------------------------------------------------------------------------
# POR QUÉ HACE FALTA ESTE CONTROL.
# La Sección 4.3 deja el test pareado rechazando la igualdad de márgenes, de
# modo que la lectura natural sería "el margen se expande hacia el cierre".
# Pero la Sección 4.2 muestra que de las cuatro temporadas solo 2021/22 es
# significativa, y que la diferencia decae hasta desaparecer en las
# siguientes. Eso abre la posibilidad de que el resultado agregado no sea un
# comportamiento del período sino un rasgo de una única temporada arrastrado
# al promedio.
#
# Son dos afirmaciones muy distintas para la tesina:
#   (a) "Bet365 amplía el margen hacia el cierre"  -> comportamiento estructural
#   (b) "En 2021/22 Bet365 ampliaba el margen hacia el cierre, práctica que
#        desaparece en las temporadas siguientes"  -> episodio acotado
#
# La única forma de decidir entre (a) y (b) es reestimar el test sin la
# temporada que concentra el efecto. Si al sacarla el resultado desaparece,
# corresponde (b). Es un control que cualquier lector atento de la tabla por
# temporada va a pedir, así que conviene traerlo resuelto.

cat("\n--- Robustez temporal: test pareado excluyendo 2021/22 ---\n")

pareado_sin_2122 <- pareado_robusto |> filter(temporada != "2021/22")

t_sin <- t.test(pareado_sin_2122$overround_c, pareado_sin_2122$overround,
                paired = TRUE)
w_sin <- wilcox.test(pareado_sin_2122$overround_c, pareado_sin_2122$overround,
                     paired = TRUE, exact = FALSE)

robustez_temporal <- tibble(
  version    = c("Sin extremos (4 temporadas)",
                 "Sin extremos y sin 2021/22"),
  n          = c(nrow(pareado_robusto), nrow(pareado_sin_2122)),
  dif_pp     = c(mean(pareado_robusto$dif), mean(pareado_sin_2122$dif)) * 100,
  t          = c(t_rob$statistic, t_sin$statistic),
  p_t        = c(t_rob$p.value,   t_sin$p.value),
  p_wilcoxon = c(w_rob$p.value,   w_sin$p.value)
)

robustez_temporal |>
  mutate(across(where(is.numeric), ~ round(., 5))) |>
  print(width = Inf)

efecto_sobrevive <- t_sin$p.value < 0.05 && w_sin$p.value < 0.05
efecto_desaparece <- t_sin$p.value >= 0.05 && w_sin$p.value >= 0.05

if (efecto_desaparece) {
  cat("\nAl excluir 2021/22 el efecto DESAPARECE en ambos tests.\n")
  cat("La expansión del margen hacia el cierre NO es un patrón del período\n")
  cat("2021-2025 sino un rasgo de la temporada 2021/22. Debe reportarse como\n")
  cat("episodio acotado y no como comportamiento estructural de la casa de\n")
  cat("apuestas. Es coherente con la serie por temporada, donde la diferencia\n")
  cat("decae de +0,119 pp hasta valores indistinguibles de cero.\n")
  cat("\nHALLAZGO REPORTABLE: no es que no pase nada, es que la política de\n")
  cat("márgenes CAMBIÓ durante el período analizado. Eso es más informativo\n")
  cat("que un promedio que mezcla ambos regímenes.\n")
} else if (efecto_sobrevive) {
  cat("\nEl efecto SE MANTIENE al excluir 2021/22 en ambos tests.\n")
  cat("La expansión del margen es un patrón del período y no un episodio,\n")
  cat("aunque su intensidad varíe entre temporadas.\n")
} else {
  cat("\nLos tests discrepan al excluir 2021/22. Reportar ambos y describir el\n")
  cat("resultado como sensible a la composición temporal de la muestra.\n")
}


# -----------------------------------------------------------------------------
# 4.5 CONCLUSIÓN DE LA COMPARACIÓN PAREADA
# -----------------------------------------------------------------------------

dif_media_pp   <- mean(pareado$dif) * 100
dif_mediana_pp <- median(pareado$dif) * 100
or_medio_pre   <- mean(pareado$overround) * 100
dif_relativa   <- dif_media_pp / or_medio_pre * 100
sig_t          <- t_pareado$p.value < 0.05
sig_w          <- w_pareado$p.value < 0.05
coinciden      <- sig_t == sig_w

cat("\n", strrep("-", 78), "\n", sep = "")
cat("CONCLUSIÓN — SECCIÓN 4\n")
cat(strrep("-", 78), "\n", sep = "")

# --- Dirección y magnitud ---
direccion <- if (dif_media_pp < 0) "se COMPRIME" else "se EXPANDE"
cat(sprintf("El margen %s hacia el cierre: %+.4f pp de media (%+.4f pp de\n",
            direccion, dif_media_pp, dif_mediana_pp))
cat(sprintf("mediana) sobre un overround pre-cierre de %.4f%%, es decir un\n",
            or_medio_pre))
cat(sprintf("cambio relativo del %+.2f%%.\n", dif_relativa))

# --- Significatividad ---
if (coinciden) {
  if (sig_t) {
    cat(sprintf("\nAmbos tests rechazan la igualdad de márgenes (t: p = %.6f; Wilcoxon: p = %.6f).\n",
                t_pareado$p.value, w_pareado$p.value))
    cat("La conclusión no depende del supuesto de normalidad.\n")
  } else {
    cat(sprintf("\nNingún test rechaza la igualdad de márgenes (t: p = %.6f; Wilcoxon: p = %.6f).\n",
                t_pareado$p.value, w_pareado$p.value))
    cat("El margen no cambia de forma detectable entre ambos momentos.\n")
  }
} else {
  cat(sprintf("\nLos tests discrepan (t: p = %.6f; Wilcoxon: p = %.6f).\n",
              t_pareado$p.value, w_pareado$p.value))
  if (hay_extremos && cambia_conclusion) {
    cat("La Sección 4.3 identifica la causa: un valor extremo del overround.\n")
    cat("Excluyéndolo, ambos tests coinciden. Se reporta la versión completa\n")
    cat("y la robusta, y se privilegia esta última en la interpretación.\n")
  } else {
    cat("La causa es la asimetría de la distribución de diferencias. En ese\n")
    cat("caso corresponde privilegiar el Wilcoxon, que no supone normalidad,\n")
    cat("y reportar ambos.\n")
  }
}

cat(sprintf("\nIC 95%% de la diferencia media: [%+.4f ; %+.4f] pp.\n",
            t_pareado$conf.int[1] * 100, t_pareado$conf.int[2] * 100))
if (t_pareado$conf.int[1] * 100 > 0 || t_pareado$conf.int[2] * 100 < 0) {
  cat("El intervalo no contiene el cero.\n")
} else {
  cat("El intervalo contiene el cero.\n")
}

# --- Lectura actuarial ---
# La decisión se toma mirando AMBOS tests, no solo el t. Si discrepan, manda
# el Wilcoxon, y entonces la dirección se lee por la MEDIANA, que es el
# estadístico que ese test efectivamente contrasta.
cat("\nLECTURA ACTUARIAL\n")

algun_test_rechaza <- sig_t || sig_w
ambos_rechazan     <- sig_t && sig_w

if (!algun_test_rechaza) {
  cat("Ningún test detecta cambio en el margen entre ambos momentos.\n")
  cat("El overround opera como recargo estructural fijo y no como componente\n")
  cat("que absorba incertidumbre variable. Es un resultado con contenido\n")
  cat("propio: acota la interpretación del margen como recargo de seguridad.\n")

} else {

  referencia <- if (ambos_rechazan) dif_media_pp else dif_mediana_pp
  etiqueta   <- if (ambos_rechazan) "media" else "mediana"

  if (!ambos_rechazan) {
    cat(sprintf("Solo el Wilcoxon rechaza. Se lee la diferencia por su %s\n", etiqueta))
    cat(sprintf("(%+.4f pp), que es el estadístico que ese test contrasta.\n\n", referencia))
  }

  # ALCANCE ANTES QUE DIRECCIÓN. La Sección 4.4 ya estableció si el resultado
  # agregado sobrevive al excluir la temporada que concentra el efecto. Si no
  # sobrevive, describir el promedio como un comportamiento del período es
  # incorrecto por más que el test sobre las cuatro temporadas rechace: el
  # promedio mezcla dos regímenes. Por eso acá se declara primero el alcance
  # temporal y recién después la dirección y el mecanismo.
  dif_2122 <- pareado_temporada$dif_pp[pareado_temporada$temporada == "2021/22"]
  if (efecto_desaparece) {
    alcance <- "en la temporada 2021/22"
  } else {
    alcance <- "en el período 2021-2025"
  }

  if (efecto_desaparece) {
    cat("ALCANCE DEL RESULTADO. La Sección 4.4 mostró que este efecto NO\n")
    cat("sobrevive al excluir 2021/22. El promedio de las cuatro temporadas\n")
    cat("mezcla dos regímenes distintos, de modo que NO corresponde atribuirle\n")
    cat("a la casa un comportamiento estable en el período.\n")
    cat(sprintf("Lo que sí se puede afirmar: en 2021/22 la diferencia fue de %+.4f pp\n",
                dif_2122))
    cat("y decae hasta valores indistinguibles de cero en las temporadas\n")
    cat("siguientes. Lo que sigue describe ESE episodio, no el período.\n\n")
  }

  if (referencia < 0) {
    cat(sprintf("El margen se COMPRIME hacia el cierre %s. Es consistente\n", alcance))
    cat("con que una porción del overround remunere incertidumbre y no\n")
    cat("únicamente margen comercial: al acercarse el partido la información\n")
    cat("relevante se hace pública y el recargo defensivo pierde justificación.\n")
    cat("Es la distinción entre recargo de seguridad y recargo de gastos y\n")
    cat("beneficio en tarificación de seguros no vida.\n")
  } else {
    cat(sprintf("El margen se EXPANDE hacia el cierre %s. NO responde a la\n", alcance))
    cat("lectura de recargo por incertidumbre: si el recargo cubriera\n")
    cat("incertidumbre, debería caer a medida que esta se resuelve.\n\n")
    cat("Interpretaciones compatibles con una expansión:\n")
    cat("  - gestión de exposición sobre el volumen ya comprometido;\n")
    cat("  - protección ante flujo de apuestas informado en las horas previas,\n")
    cat("    que es el mecanismo de Shin (1993) operando en sentido inverso al\n")
    cat("    esperado: el riesgo de agentes informados sería MAYOR, no menor,\n")
    cat("    cerca del cierre.\n\n")
    if (efecto_desaparece) {
      cat("REPORTAR COMO CAMBIO DE RÉGIMEN, NO COMO RASGO DEL PERÍODO. Estos\n")
      cat("mecanismos explicarían por qué la casa aplicaba esa política en\n")
      cat("2021/22; el hallazgo reportable es que la abandonó después.\n")
    } else {
      cat("REPORTAR COMO HALLAZGO Y NO FORZAR LA INTERPRETACIÓN.\n")
    }
    cat(sprintf("La magnitud es pequeña: %+.4f pp sobre un margen de %.2f%%, es\n",
                referencia, or_medio_pre))
    cat(sprintf("decir %+.2f%% en términos relativos. Es un cambio detectable pero\n",
                referencia / or_medio_pre * 100))
    cat("económicamente marginal, y conviene decirlo así.\n")
  }
}

# --- Estabilidad entre temporadas ---
signos <- sign(pareado_temporada$dif_pp)
cat("\nESTABILIDAD ENTRE TEMPORADAS\n")
if (length(unique(signos)) == 1) {
  cat("El signo de la diferencia se mantiene en las cuatro temporadas:\n")
  cat("el patrón es estable y no lo produce una temporada atípica.\n")
} else {
  cat("El signo de la diferencia NO es constante entre temporadas.\n")
  cat("El promedio agregado esconde comportamientos distintos: reportar la\n")
  cat("desagregación y no solo el dato global. Verificar en la tabla si el\n")
  cat("efecto global lo produce una sola temporada.\n")
}
print(pareado_temporada)

temporadas_sig <- sum(pareado_temporada$p_valor < 0.05)
cat(sprintf("\nTemporadas con diferencia significativa: %d de %d.\n",
            temporadas_sig, nrow(pareado_temporada)))

# La Sección 4.4 ya resolvió empíricamente si el efecto sobrevive sin la
# temporada dominante. Acá se traduce ese resultado a la redacción.
if (temporadas_sig == 1 && efecto_desaparece) {
  cat("Una sola temporada concentra el efecto y, segun la Seccion 4.4, el\n")
  cat("resultado NO sobrevive al excluirla.\n")
  cat("\nREDACCION CORRECTA: no afirmar que el margen se expande hacia el\n")
  cat("cierre. Afirmar que la casa aplicaba esa politica en 2021/22 y la\n")
  cat("abandono despues. El hallazgo es el CAMBIO DE REGIMEN, no el promedio.\n")
} else if (temporadas_sig == 1 && efecto_sobrevive) {
  cat("Una sola temporada alcanza significacion por separado, pero el efecto\n")
  cat("agregado sobrevive al excluirla (Seccion 4.4): el patron existe en\n")
  cat("todo el periodo y 2021/22 solo es donde se ve con mas claridad.\n")
} else {
  cat("Revisar la Seccion 4.4 para decidir si el resultado agregado describe\n")
  cat("el periodo completo o esta dominado por una temporada.\n")
}



# -----------------------------------------------------------------------------
# 5. GRÁFICOS DEL OVERROUND
# -----------------------------------------------------------------------------

# --- Gráfico 1: distribución del overround en ambas series ---
g_overround_dist <- datos_mercado |>
  select(overround, overround_c) |>
  pivot_longer(everything(), names_to = "serie", values_to = "valor") |>
  filter(!is.na(valor)) |>
  mutate(serie = if_else(serie == "overround",
                         "Pre-cierre (B365H/D/A)",
                         "Cierre (B365CH/CD/CA)")) |>
  ggplot(aes(x = valor * 100, fill = serie)) +
  geom_histogram(binwidth = 0.2, color = "white", alpha = 0.75,
                 position = "identity") +
  facet_wrap(~serie, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Pre-cierre (B365H/D/A)"   = "#2c3e50",
                               "Cierre (B365CH/CD/CA)"    = "#16a085")) +
  labs(
    title    = "Distribución del overround — Bet365",
    subtitle = sprintf("Premier League 2021/22–2024/25 | N = %d partidos",
                       nrow(datos_mercado)),
    x        = "Overround (%)",
    y        = "Frecuencia",
    fill     = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(g_overround_dist)

# --- Gráfico 2: evolución del overround medio por temporada ---
g_overround_temp <- overround_temporada |>
  select(temporada, pre_cierre, cierre) |>
  pivot_longer(cols = c(pre_cierre, cierre),
               names_to = "serie", values_to = "overround") |>
  mutate(serie = if_else(serie == "pre_cierre", "Pre-cierre", "Cierre")) |>
  ggplot(aes(x = temporada, y = overround, color = serie, group = serie)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Pre-cierre" = "#2c3e50",
                                "Cierre"     = "#16a085")) +
  labs(
    title    = "Evolución del overround medio por temporada",
    subtitle = "Bet365 — Premier League",
    x        = NULL,
    y        = "Overround medio (%)",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g_overround_temp)

# --- Gráfico 3: overround según la probabilidad del favorito ---
# ¿La casa aplica el mismo margen a favoritos y a longshots?
g_overround_fav <- datos_mercado |>
  mutate(
    pi_max      = pmax(pi_H, pi_A),   # se ignora el empate: distorsiona el eje
    bin_mercado = ntile(pi_max, 10)
  ) |>
  summarise(
    pi_max_med       = mean(pi_max),
    overround_med    = mean(overround,   na.rm = TRUE) * 100,
    overround_c_med  = mean(overround_c, na.rm = TRUE) * 100,
    .by = bin_mercado
  ) |>
  pivot_longer(cols = c(overround_med, overround_c_med),
               names_to = "serie", values_to = "overround") |>
  mutate(serie = if_else(serie == "overround_med", "Pre-cierre", "Cierre")) |>
  ggplot(aes(x = pi_max_med * 100, y = overround, color = serie)) +
  geom_point(size = 3.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_color_manual(values = c("Pre-cierre" = "#2c3e50",
                                "Cierre"     = "#16a085")) +
  labs(
    title    = "Overround según la probabilidad del favorito",
    subtitle = "Cada punto es un decil de partidos | Premier League 2021/22–2024/25",
    x        = "Probabilidad implícita del favorito (promedio del decil, %)",
    y        = "Overround promedio (%)",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g_overround_fav)


# -----------------------------------------------------------------------------
# 6. FORMATO LARGO PARA EL ANÁLISIS DE MERCADO
# -----------------------------------------------------------------------------
# Se apilan los tres resultados posibles (H, D, A) de cada partido en filas
# separadas: de ≈1.520 partidos se pasa a ≈4.560 filas resultado × partido.
# Para cada fila queda registrada la probabilidad implícita normalizada de
# ambas series, la cuota decimal correspondiente, y el indicador binario
# de ocurrencia.
#
# ATENCIÓN: estas ≈4.560 filas NO son independientes entre sí (las tres de un
# mismo partido suman 1). Por eso el objeto largo se usa solo para armar
# gráficos y tablas descriptivas; la INFERENCIA se hace en la Sección 7,
# corriendo una regresión por categoría, donde cada regresión tiene
# exactamente una observación por partido.

mercado_largo <- datos_mercado |>
  select(temporada, date, home_team, away_team, ftr,
         pi_H, pi_D, pi_A,
         piC_H, piC_D, piC_A,
         cuota_H  = b365h,  cuota_D  = b365d,  cuota_A  = b365a,
         cuotaC_H = b365ch, cuotaC_D = b365cd, cuotaC_A = b365ca) |>
  pivot_longer(
    cols      = c(pi_H, pi_D, pi_A, piC_H, piC_D, piC_A,
                  cuota_H, cuota_D, cuota_A, cuotaC_H, cuotaC_D, cuotaC_A),
    names_to  = c(".value", "resultado"),
    names_sep = "_"
  ) |>
  rename(
    pi_precierre    = pi,
    pi_cierre       = piC,
    cuota_precierre = cuota,
    cuota_cierre    = cuotaC
  ) |>
  mutate(ocurrio = as.integer(ftr == resultado))

cat("\n=== FORMATO LARGO (análisis de mercado) ===\n")
cat("Filas:", nrow(mercado_largo), "=", nrow(datos_mercado), "partidos × 3 resultados\n")
cat("\nDistribución de resultados observados:\n")
mercado_largo |>
  summarise(n = n(), ocurrencias = sum(ocurrio),
            frecuencia = round(mean(ocurrio), 4), .by = resultado) |>
  print()



# -----------------------------------------------------------------------------
# 7. TEST FORMAL DE FLB: REGRESIÓN LOGÍSTICA DE CALIBRACIÓN
# -----------------------------------------------------------------------------
# 7.1 REGRESIONES POR CATEGORÍA Y POR SERIE
# -----------------------------------------------------------------------------
# Especificación, por cada categoría (H, D, A) y cada serie, POR SEPARADO:
#
#     logit( P(ocurrio = 1) ) = a + b · logit(π)
#
# donde π es la probabilidad implícita NORMALIZADA de esa categoría.
#
# Bajo eficiencia perfecta:  a = 0  y  b = 1.
#
# El FLB implica b > 1: las probabilidades verdaderas son más extremas que las
# implícitas —los favoritos ganan MÁS de lo que sugiere su cuota y los
# longshots MENOS—, de modo que la relación entre logit verdadero y logit
# implícito tiene pendiente mayor que uno.
#
# Contraste relevante: H0: b = 1, mediante el estadístico de Wald
#     z = (b − 1) / se(b)
# NO el p-valor por defecto de summary(), que contrasta H0: b = 0 y en este
# contexto no dice nada (b = 0 significaría que la cuota no informa nada,
# hipótesis que nadie sostiene).
#
# Se agrega además el contraste conjunto H0: (a, b) = (0, 1) mediante el
# estadístico de Wald con 2 grados de libertad, que es el test completo de
# calibración: rechaza si la cuota está sesgada por nivel, por pendiente, o
# por ambos.
#
# POR QUÉ TRES REGRESIONES SEPARADAS Y NO UNA SOLA:
# las tres categorías de un mismo partido no son independientes (suman 1).
# Corriendo una regresión por categoría, cada una tiene exactamente una
# observación por partido y constituye una muestra genuinamente independiente
# de ≈1.520 casos. Además es informativo por sí mismo: no hay razón a priori
# para que el sesgo opere igual sobre locales, empates y visitantes.

regresion_calibracion <- function(df, categoria, columna_pi, etiqueta_serie) {

  d <- df |>
    filter(resultado == categoria) |>
    transmute(ocurrio, p = .data[[columna_pi]]) |>
    filter(!is.na(p), p > 0, p < 1) |>
    mutate(logit_pi = log(p / (1 - p)))

  modelo_glm <- glm(ocurrio ~ logit_pi,
                    family = binomial(link = "logit"),
                    data   = d)

  a  <- unname(coef(modelo_glm)[1])
  b  <- unname(coef(modelo_glm)[2])
  V  <- vcov(modelo_glm)
  se_a <- sqrt(V[1, 1])
  se_b <- sqrt(V[2, 2])

  # Wald sobre H0: b = 1
  z_b <- (b - 1) / se_b
  p_b <- 2 * pnorm(-abs(z_b))

  # Wald conjunto sobre H0: (a, b) = (0, 1)
  theta <- c(a, b - 1)
  W     <- as.numeric(t(theta) %*% solve(V) %*% theta)
  p_W   <- pchisq(W, df = 2, lower.tail = FALSE)

  tibble(
    serie          = etiqueta_serie,
    categoria      = categoria,
    n              = nrow(d),
    intercepto     = a,
    se_intercepto  = se_a,
    pendiente      = b,
    se_pendiente   = se_b,
    ic95_inf       = b - qnorm(0.975) * se_b,
    ic95_sup       = b + qnorm(0.975) * se_b,
    z_pendiente    = z_b,      # contra H0: b = 1
    p_pendiente    = p_b,
    wald_conjunto  = W,        # contra H0: (a, b) = (0, 1)
    p_conjunto     = p_W
  )
}

categorias <- c("H", "D", "A")

flb_regresiones <- bind_rows(
  map_dfr(categorias, ~ regresion_calibracion(mercado_largo, .x,
                                              "pi_precierre", "Pre-cierre")),
  map_dfr(categorias, ~ regresion_calibracion(mercado_largo, .x,
                                              "pi_cierre",    "Cierre"))
) |>
  mutate(
    etiqueta_cat = case_when(categoria == "H" ~ "Local (H)",
                             categoria == "D" ~ "Empate (D)",
                             categoria == "A" ~ "Visitante (A)"),
    conclusion = case_when(
      p_pendiente >= 0.05                    ~ "No se rechaza b = 1",
      p_pendiente <  0.05 & pendiente > 1    ~ "b > 1: consistente con FLB",
      TRUE                                   ~ "b < 1: sesgo inverso al FLB"
    )
  )

cat("\n", strrep("=", 78), "\n", sep = "")
cat("REGRESIÓN LOGÍSTICA DE CALIBRACIÓN — TEST FORMAL DE FLB\n")
cat(strrep("=", 78), "\n", sep = "")
cat("Especificación: logit(P(ocurrio)) = a + b · logit(pi)\n")
cat("Eficiencia perfecta: a = 0, b = 1  |  FLB: b > 1\n")
cat("Muestra: 4 temporadas, una observación por partido en cada regresión\n\n")

flb_regresiones |>
  select(serie, etiqueta_cat, n, pendiente, se_pendiente,
         ic95_inf, ic95_sup, z_pendiente, p_pendiente) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

cat("\n--- Intercepto y contraste conjunto H0: (a, b) = (0, 1) ---\n")
flb_regresiones |>
  select(serie, etiqueta_cat, intercepto, se_intercepto,
         wald_conjunto, p_conjunto) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

cat("\n--- Conclusión por regresión ---\n")
flb_regresiones |>
  select(serie, etiqueta_cat, pendiente, p_pendiente, conclusion) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

# --- Gráfico 4: pendientes estimadas con intervalo de confianza ---
g_flb_pendientes <- flb_regresiones |>
  mutate(etiqueta_cat = factor(etiqueta_cat,
                               levels = c("Local (H)", "Empate (D)", "Visitante (A)"))) |>
  ggplot(aes(x = etiqueta_cat, y = pendiente, color = serie)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup),
                  position = position_dodge(width = 0.4), size = 0.7) +
  scale_color_manual(values = c("Pre-cierre" = "#2c3e50",
                                "Cierre"     = "#16a085")) +
  labs(
    title    = "Pendiente de la regresión de calibración con IC del 95%",
    subtitle = "Línea de referencia b = 1 (eficiencia) | b > 1 indica favourite-longshot bias",
    x        = NULL,
    y        = "Pendiente estimada (b)",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g_flb_pendientes)



# -----------------------------------------------------------------------------
# 7.2 TEST AGREGADO CON ERRORES ESTÁNDAR AGRUPADOS POR PARTIDO
# -----------------------------------------------------------------------------
# POR QUÉ HACE FALTA. Las seis regresiones de la Sección 7.1 dan seis
# respuestas parciales y correlacionadas entre sí: las dos series comparten
# los mismos partidos, y dentro de cada partido las categorías no son
# independientes. Con p-valores individuales no hay forma limpia de decidir
# si el conjunto de la evidencia apunta o no a un sesgo.
#
# LA SOLUCIÓN ESTÁNDAR. Apilar las tres categorías en una sola regresión con
# interceptos propios por categoría y una PENDIENTE COMÚN, y corregir los
# errores estándar por el agrupamiento a nivel partido. El estimador de
# varianza agrupado (cluster-robust) no exige independencia dentro del
# partido: solo entre partidos, que es exactamente la estructura del dato.
#
# Especificación:
#   logit(P(ocurrio)) = a_H·1{H} + a_D·1{D} + a_A·1{A} + b·logit(pi)
# Contraste: H0: b = 1
#
# Los interceptos por categoría absorben las diferencias de frecuencia base
# entre local, empate y visitante, de modo que b mide exclusivamente la
# distorsión de calibración común a las tres.
#
# Resultado: un único contraste con inferencia correcta, que es el número que
# encabeza el capítulo. Las seis regresiones quedan como desagregación.

# --- Estimador de varianza agrupado, implementado sin paquetes externos ---
# V_CR = (X'WX)^-1 [ suma_g U_g U_g' ] (X'WX)^-1
# con U_g = suma sobre i en g de x_i (y_i - mu_i): el score agregado del
# partido g.
vcov_agrupado <- function(modelo, grupo) {
  X <- model.matrix(modelo)
  e <- residuals(modelo, type = "response")     # y - mu
  U <- rowsum(X * e, factor(grupo))             # score por partido
  pan   <- vcov(modelo)                         # (X'WX)^-1
  carne <- crossprod(U)
  G <- nlevels(factor(grupo)); n <- nrow(X); k <- ncol(X)
  correccion <- (G / (G - 1)) * ((n - 1) / (n - k))
  pan %*% carne %*% pan * correccion
}

test_agregado <- function(df, columna_pi, etiqueta_serie) {

  d <- df |>
    mutate(id_partido = paste(date, home_team, away_team)) |>
    transmute(id_partido, ocurrio,
              categoria = factor(resultado, levels = c("H", "D", "A")),
              p = .data[[columna_pi]]) |>
    filter(!is.na(p), p > 0, p < 1) |>
    mutate(logit_pi = log(p / (1 - p)))

  m <- glm(ocurrio ~ 0 + categoria + logit_pi,
           family = binomial(link = "logit"), data = d)

  V_cr <- vcov_agrupado(m, d$id_partido)

  b      <- unname(coef(m)["logit_pi"])
  se_ing <- sqrt(vcov(m)["logit_pi", "logit_pi"])   # ingenuo (incorrecto acá)
  se_cr  <- sqrt(V_cr["logit_pi", "logit_pi"])      # agrupado (correcto)
  z      <- (b - 1) / se_cr

  tibble(
    serie         = etiqueta_serie,
    partidos      = n_distinct(d$id_partido),
    observaciones = nrow(d),
    pendiente     = b,
    se_ingenuo    = se_ing,
    se_agrupado   = se_cr,
    ic95_inf      = b - qnorm(0.975) * se_cr,
    ic95_sup      = b + qnorm(0.975) * se_cr,
    z             = z,
    p_valor       = 2 * pnorm(-abs(z))
  )
}

flb_agregado <- bind_rows(
  test_agregado(mercado_largo, "pi_precierre", "Pre-cierre"),
  test_agregado(mercado_largo, "pi_cierre",    "Cierre")
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("TEST AGREGADO DE FLB — PENDIENTE COMÚN, EE AGRUPADOS POR PARTIDO\n")
cat(strrep("=", 78), "\n", sep = "")
cat("Interceptos por categoría + pendiente común | H0: b = 1\n")
cat("Varianza agrupada a nivel partido\n\n")

flb_agregado |>
  select(serie, partidos, observaciones, pendiente,
         se_ingenuo, se_agrupado, ic95_inf, ic95_sup, z, p_valor) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(width = Inf)

cat("\n--- Efecto de ignorar el agrupamiento ---\n")
cat("se_ingenuo es el error estándar que saldría de tratar las 4.560 filas\n")
cat("como independientes. El cociente muestra cuánto se distorsionaría la\n")
cat("inferencia al no corregir:\n")
flb_agregado |>
  transmute(serie, cociente = round(se_agrupado / se_ingenuo, 3)) |>
  print(width = Inf)
cat("(cociente > 1 = el EE ingenuo subestima la incertidumbre)\n")

for (i in seq_len(nrow(flb_agregado))) {
  f <- flb_agregado[i, ]
  cat(sprintf("\n%s: b = %.4f | IC 95%% [%.4f ; %.4f] | p = %.4f\n",
              f$serie, f$pendiente, f$ic95_inf, f$ic95_sup, f$p_valor))
  if (f$p_valor < 0.05 && f$pendiente > 1) {
    cat("  Se rechaza b = 1 en dirección consistente con FLB.\n")
  } else if (f$p_valor < 0.05) {
    cat("  Se rechaza b = 1 en dirección inversa al FLB.\n")
  } else {
    cat("  No se rechaza b = 1.\n")
  }
  cat("  El intervalo es compatible con ausencia de sesgo (b = 1) y\n")
  cat(sprintf("  también con un sesgo de hasta b = %.2f. El estudio no los\n", f$ic95_sup))
  cat("  distingue: es la afirmación precisa que corresponde reportar.\n")
}



# -----------------------------------------------------------------------------
# 7.3 CONCLUSIÓN DEL TEST DE FLB
# -----------------------------------------------------------------------------

cat("\n", strrep("-", 78), "\n", sep = "")
cat("CONCLUSIÓN — SECCIÓN 7\n")
cat(strrep("-", 78), "\n", sep = "")

# --- 1. Resultado principal: el test agregado ---
cat("\nRESULTADO PRINCIPAL (test agregado, Sección 7.2)\n")
for (i in seq_len(nrow(flb_agregado))) {
  f <- flb_agregado[i, ]
  cat(sprintf("  %-11s b = %.4f  IC 95%% [%.3f ; %.3f]  p = %.4f  -> %s\n",
              f$serie, f$pendiente, f$ic95_inf, f$ic95_sup, f$p_valor,
              if (f$p_valor < 0.05 && f$pendiente > 1) "rechaza, FLB"
              else if (f$p_valor < 0.05)               "rechaza, inverso"
              else                                     "no rechaza"))
}

# --- 2. Desagregación: las seis regresiones ---
sintesis_flb <- flb_regresiones |>
  summarise(
    regresiones = n(),
    rechazan    = sum(p_pendiente < 0.05),
    con_flb     = sum(p_pendiente < 0.05 & pendiente > 1),
    inverso     = sum(p_pendiente < 0.05 & pendiente < 1),
    .by = serie
  )

cat("\n--- Desagregación: rechazos de H0: b = 1 por serie (alfa = 0,05) ---\n")
print(sintesis_flb, width = Inf)

total_rechazos <- sum(sintesis_flb$rechazan)
total_flb      <- sum(sintesis_flb$con_flb)
total_inverso  <- sum(sintesis_flb$inverso)

if (total_rechazos == 0) {
  cat("\nNinguna regresión individual rechaza la eficiencia.\n")
} else {
  cat(sprintf("\nSe rechaza H0: b = 1 en %d de las 6 regresiones", total_rechazos))
  cat(sprintf(" (%d con b > 1, %d con b < 1).\n", total_flb, total_inverso))
  if (total_rechazos <= 2) {
    cat("CAUTELA: con seis contrastes al 5%, uno o dos rechazos aislados son\n")
    cat("compatibles con el azar. Pero los seis contrastes NO son independientes\n")
    cat("—comparten partidos y categorías—, de modo que la corrección por\n")
    cat("comparaciones múltiples tampoco corresponde aplicarla como si lo\n")
    cat("fueran. Esa es precisamente la razón por la que el test agregado de\n")
    cat("la Sección 7.2 es el instrumento de referencia.\n")
  }
}

# --- 3. Síntesis direccional ---
# La dispersión de logit(pi) determina qué categorías pueden informar algo.
# La Sección 7.4 la cuantifica en términos de potencia; acá se usa solo para
# separar las categorías con y sin capacidad de estimar una pendiente.
dispersion_cat <- mercado_largo |>
  filter(pi_precierre > 0, pi_precierre < 1) |>
  summarise(sd_logit = sd(log(pi_precierre / (1 - pi_precierre))),
            .by = resultado)

cat("\n--- Síntesis direccional: signo de las desviaciones respecto de b = 1 ---\n")
flb_regresiones |>
  left_join(dispersion_cat, by = c("categoria" = "resultado")) |>
  transmute(serie, etiqueta_cat,
            sd_logit  = round(sd_logit, 3),
            pendiente = round(pendiente, 4),
            supera_1  = pendiente > 1) |>
  print(n = Inf, width = Inf)

n_supera_total <- sum(flb_regresiones$pendiente > 1)
reg_HA         <- flb_regresiones |> filter(categoria != "D")
n_supera_HA    <- sum(reg_HA$pendiente > 1)

cat(sprintf("\nPendientes por encima de 1: %d de %d en total; %d de %d excluyendo\n",
            n_supera_total, nrow(flb_regresiones), n_supera_HA, nrow(reg_HA)))
cat("los empates, cuya dispersión es aproximadamente un tercio de la de las\n")
cat("otras dos categorías.\n")

if (n_supera_HA == nrow(reg_HA)) {
  cat("\nTODAS las estimaciones con dispersión suficiente exceden 1. La\n")
  cat("consistencia direccional es un dato a reportar. Ahora bien, NO se\n")
  cat("acumula en significación: las estimaciones no son independientes\n")
  cat("entre sí, de modo que no corresponde aplicarles un test de signos ni\n")
  cat("tratarlas como confirmaciones separadas. Lo que sintetiza esa\n")
  cat("evidencia con inferencia válida es el test agregado de la 7.2.\n")
}

# --- 4. Concordancia entre series ---
concordancia <- flb_regresiones |>
  mutate(rechaza = p_pendiente < 0.05) |>
  select(etiqueta_cat, serie, rechaza, pendiente) |>
  pivot_wider(names_from = serie, values_from = c(rechaza, pendiente)) |>
  mutate(
    coinciden     = `rechaza_Pre-cierre` == rechaza_Cierre,
    dif_pendiente = round(pendiente_Cierre - `pendiente_Pre-cierre`, 4)
  )

cat("\n--- Concordancia entre series por categoría ---\n")
print(concordancia, width = Inf)

if (all(concordancia$coinciden)) {
  cat("\nLas dos series coinciden en las tres categorías: la conclusión sobre\n")
  cat("el sesgo no depende de qué serie de cuotas se utilice.\n")
} else {
  cat("\nLas series no coinciden en todas las categorías EN TÉRMINOS DE\n")
  cat("SIGNIFICACIÓN. Antes de leerlo como contradicción, mirar la columna\n")
  cat("dif_pendiente: si las pendientes son casi iguales, lo que difiere es\n")
  cat("de qué lado de 0,05 cae el p-valor, no la estimación. Reportar ambas\n")
  cat("series y describirlo como diferencia de precisión, no de resultado.\n")
  cat("Lo que NO corresponde es elegir la serie más conveniente.\n")
}

# --- 5. Contraste conjunto (a, b) = (0, 1) ---
cat("\n--- Contraste conjunto de calibración H0: (a, b) = (0, 1) ---\n")
n_conjunto <- sum(flb_regresiones$p_conjunto < 0.05)
cat(sprintf("Rechaza en %d de 6 regresiones.\n", n_conjunto))
if (n_conjunto > total_rechazos) {
  cat("El contraste conjunto rechaza más veces que el de la pendiente sola.\n")
  cat("Indica desajuste de NIVEL (intercepto distinto de cero) con pendiente\n")
  cat("compatible con 1: el mercado tendría un sesgo constante, no del tipo\n")
  cat("favourite-longshot. Es un hallazgo distinto y hay que nombrarlo así.\n")
} else if (n_conjunto == 0) {
  cat("Ningún rechazo del contraste conjunto: no hay evidencia de desajuste\n")
  cat("de nivel. Los interceptos son compatibles con cero, de modo que el\n")
  cat("mercado no muestra un sesgo constante superpuesto al de pendiente.\n")
}



# -----------------------------------------------------------------------------
# 7.4 POTENCIA DEL TEST (SIMULACIÓN SOBRE LOS DATOS OBSERVADOS)
# -----------------------------------------------------------------------------
# 7.4.1 SIMULACIÓN
# -----------------------------------------------------------------------------
# El objetivo declarado del rediseño era pasar de "no se detectó FLB" a
# "no se detectó FLB con un diseño con potencia razonable para detectarlo".
# Esa segunda afirmación es más fuerte, pero solo es legítima si la potencia
# se calcula: si no, es una promesa.
#
# Procedimiento (simulación paramétrica sobre la distribución empírica de π):
#   1. Se toma el vector observado de probabilidades implícitas de la categoría.
#   2. Se genera un resultado artificial bajo la hipótesis alternativa
#      P(ocurre) = logit⁻¹( b_alt · logit(π) ), con b_alt > 1 (FLB de tamaño b_alt).
#   3. Se estima la regresión y se aplica el test de Wald contra H0: b = 1.
#   4. La potencia es la proporción de réplicas en que se rechaza.
#
# Se usa la distribución de π de la serie pre-cierre. La de cierre es casi
# idéntica en dispersión, de modo que la potencia es equivalente.
#
# ADVERTENCIA METODOLÓGICA IMPORTANTE (va al texto):
# la potencia de este test depende críticamente de la DISPERSIÓN de logit(π)
# dentro de la categoría, no solo del tamaño de muestra. El mercado cotiza
# los empates en un rango muy angosto (típicamente entre 0,15 y 0,30), con
# lo cual la varianza de logit(π) en la categoría D es una fracción de la de
# H y A. La consecuencia es estructural: en los empates el test carece de
# potencia por construcción, y un resultado no significativo en esa categoría
# NO es evidencia de ausencia de sesgo. Esto debe declararse explícitamente
# al reportar los resultados y en el capítulo de limitaciones.

potencia_flb <- function(df, categoria, columna_pi, b_alt,
                         n_sim = 500, alpha = 0.05, semilla = 904158) {

  set.seed(semilla)

  p <- df |>
    filter(resultado == categoria) |>
    pull(.data[[columna_pi]])
  p <- p[!is.na(p) & p > 0 & p < 1]

  logit_pi <- log(p / (1 - p))
  p_alt    <- 1 / (1 + exp(-(b_alt * logit_pi)))

  rechazos <- 0L
  for (s in seq_len(n_sim)) {
    y <- rbinom(length(p_alt), size = 1, prob = p_alt)
    # Si la réplica sale degenerada (todo 0 o todo 1) se descarta
    if (length(unique(y)) < 2) next
    ajuste <- suppressWarnings(
      glm(y ~ logit_pi, family = binomial(link = "logit"))
    )
    b  <- unname(coef(ajuste)[2])
    se <- sqrt(vcov(ajuste)[2, 2])
    if (is.finite(se) && se > 0) {
      if (2 * pnorm(-abs((b - 1) / se)) < alpha) rechazos <- rechazos + 1L
    }
  }

  tibble(
    categoria  = categoria,
    n          = length(p),
    sd_logit   = sd(logit_pi),
    b_alt      = b_alt,
    potencia   = rechazos / n_sim
  )
}

# Grilla de tamaños de sesgo bajo la alternativa. Se incluyen 1,15 y 1,25
# porque con una grilla más gruesa el efecto mínimo detectable de H y A
# aparecería en 1,30 cuando la potencia en 1,20 ya roza el 80%: reportar
# "solo detecta b >= 1,30" subestimaría lo que el diseño efectivamente puede.
valores_b <- c(1.10, 1.15, 1.20, 1.25, 1.30, 1.50)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("POTENCIA DEL TEST DE FLB (simulación sobre la distribución observada)\n")
cat(strrep("=", 78), "\n", sep = "")
cat("Puede tardar un minuto...\n")

tabla_potencia <- expand_grid(categoria = categorias, b_alt = valores_b) |>
  pmap_dfr(function(categoria, b_alt) {
    potencia_flb(mercado_largo, categoria, "pi_precierre", b_alt)
  }) |>
  mutate(etiqueta_cat = case_when(categoria == "H" ~ "Local (H)",
                                  categoria == "D" ~ "Empate (D)",
                                  categoria == "A" ~ "Visitante (A)"))

cat("\n--- Dispersión de logit(pi) por categoría ---\n")
cat("(a mayor dispersión, mayor potencia: es el determinante principal)\n")
tabla_potencia |>
  distinct(etiqueta_cat, n, sd_logit) |>
  mutate(sd_logit = round(sd_logit, 4)) |>
  print()

cat("\n--- Potencia estimada frente a distintos tamaños de sesgo ---\n")
tabla_potencia |>
  select(etiqueta_cat, b_alt, potencia) |>
  pivot_wider(names_from = b_alt, values_from = potencia,
              names_prefix = "b=") |>
  print(width = Inf)

cat("\nLectura: la potencia es la probabilidad de detectar un FLB de ese tamaño\n")
cat("         si efectivamente existiera. Valores bajos significan que un\n")
cat("         resultado no significativo es poco informativo.\n")

# --- Gráfico 4b: curva de potencia ---
g_potencia <- tabla_potencia |>
  mutate(etiqueta_cat = factor(etiqueta_cat,
                               levels = c("Local (H)", "Empate (D)", "Visitante (A)"))) |>
  ggplot(aes(x = b_alt, y = potencia, color = etiqueta_cat, group = etiqueta_cat)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_color_manual(values = c("Local (H)"     = "#2980b9",
                                "Empate (D)"    = "#e67e22",
                                "Visitante (A)" = "#27ae60")) +
  labs(
    title    = "Potencia del test de FLB por categoría",
    subtitle = sprintf("Simulación sobre la distribución observada de π | %d partidos | α = 0,05\nLínea discontinua: potencia del 80%%",
                       nrow(datos_mercado)),
    x        = "Tamaño del sesgo bajo la alternativa (pendiente verdadera b)",
    y        = "Potencia",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g_potencia)



# -----------------------------------------------------------------------------
# 7.4.2 CONCLUSIÓN SOBRE LA POTENCIA Y LECTURA CONJUNTA
# -----------------------------------------------------------------------------
# Este bloque convierte "no se detectó FLB" en una afirmación con contenido.
# Un no rechazo con potencia alta es evidencia de ausencia de sesgo relevante;
# un no rechazo con potencia baja no dice nada.

cat("\n", strrep("-", 78), "\n", sep = "")
cat("CONCLUSIÓN — SECCIÓN 7.4\n")
cat(strrep("-", 78), "\n", sep = "")

umbral <- 0.80

mde <- tabla_potencia |>
  filter(potencia >= umbral) |>
  slice_min(b_alt, n = 1, by = categoria) |>
  select(categoria, etiqueta_cat, b_min_detectable = b_alt, potencia)

sin_mde <- tabla_potencia |>
  filter(!categoria %in% mde$categoria) |>
  slice_max(potencia, n = 1, by = categoria) |>
  select(categoria, etiqueta_cat, b_max_probado = b_alt, potencia_max = potencia)

cat("\n--- Efecto mínimo detectable con potencia del 80% ---\n")
if (nrow(mde) > 0) {
  print(mde, width = Inf)
  cat("En esas categorías el diseño detecta sesgos de ese tamaño o mayores en\n")
  cat("al menos 8 de cada 10 muestras.\n")
}
if (nrow(sin_mde) > 0) {
  cat("\n--- Categorías que NO alcanzan el 80% en ningún valor probado ---\n")
  print(sin_mde, width = Inf)
  cat("En estas categorías el test carece de potencia adecuada incluso frente\n")
  cat("a sesgos grandes. Un no rechazo NO es evidencia de ausencia de sesgo.\n")
}

cat("\n--- Dispersión de logit(pi): el determinante de la potencia ---\n")
tabla_potencia |>
  distinct(etiqueta_cat, n, sd_logit) |>
  arrange(desc(sd_logit)) |>
  mutate(sd_logit = round(sd_logit, 4)) |>
  print(width = Inf)

cat("\nLa regresión estima una pendiente, y para estimar una pendiente hace\n")
cat("falta que el predictor varíe. Donde el mercado cotiza en un rango\n")
cat("angosto no hay palanca estadística y la potencia colapsa, por más\n")
cat("observaciones que se acumulen. Es una limitación ESTRUCTURAL del\n")
cat("mercado de empates, no un problema de tamaño de muestra.\n")

# --- LECTURA POR PRECISIÓN DEL INTERVALO ---
# Preferible a evaluar la potencia en el valor estimado. La potencia
# observada es una función monótona del p-valor y no aporta información
# adicional; el intervalo, en cambio, dice directamente qué valores de b
# el estudio descarta y cuáles no.
cat("\n--- Qué descarta y qué no descarta cada estimación ---\n")
cat("(la tabla de potencia de arriba es legítima porque evalúa valores\n")
cat(" PRE-ESPECIFICADOS de b; evaluar la potencia en el valor estimado sería\n")
cat(" una práctica desaconsejada, y el intervalo la reemplaza con ventaja)\n\n")

precision <- bind_rows(
  flb_regresiones |>
    transmute(fuente = paste(serie, "-", etiqueta_cat),
              b = pendiente, ic95_inf, ic95_sup),
  flb_agregado |>
    transmute(fuente = paste(serie, "- AGREGADO"),
              b = pendiente, ic95_inf, ic95_sup)
) |>
  mutate(
    descarta_eficiencia = ic95_inf > 1 | ic95_sup < 1,
    descarta_b_120      = ic95_sup < 1.20,
    descarta_b_130      = ic95_sup < 1.30,
    ic                  = sprintf("[%.3f ; %.3f]", ic95_inf, ic95_sup)
  ) |>
  select(fuente, b, ic, descarta_eficiencia, descarta_b_120, descarta_b_130)

precision |>
  mutate(b = round(b, 4)) |>
  print(n = Inf, width = Inf)

cat("\nLECTURA. descarta_eficiencia = TRUE significa que el estudio rechaza\n")
cat("b = 1. descarta_b_130 = TRUE significa que el estudio descarta un sesgo\n")
cat("de esa magnitud. Cuando ambas son FALSE, el intervalo es compatible\n")
cat("tanto con ausencia de sesgo como con un sesgo sustancial: el estudio\n")
cat("no los distingue, y eso es lo que hay que escribir.\n")

# --- LECTURA CONJUNTA: resultado del test cruzado con potencia ---
diagnostico <- flb_regresiones |>
  select(serie, categoria, etiqueta_cat, pendiente, p_pendiente) |>
  left_join(
    tabla_potencia |>
      filter(b_alt == 1.30) |>
      select(categoria, potencia_130 = potencia),
    by = "categoria"
  ) |>
  mutate(
    rechaza = p_pendiente < 0.05,
    lectura = case_when(
      rechaza & pendiente > 1           ~ "Sesgo detectado, consistente con FLB",
      rechaza & pendiente < 1           ~ "Sesgo detectado, de signo inverso al FLB",
      !rechaza & potencia_130 >= umbral ~ "No rechazo INFORMATIVO: descarta FLB grande",
      !rechaza & potencia_130 >= 0.50   ~ "No rechazo de fuerza intermedia",
      TRUE                              ~ "No rechazo NO CONCLUYENTE: potencia insuficiente"
    )
  )

cat("\n--- LECTURA CONJUNTA (resultado del test x potencia frente a b = 1,30) ---\n")
diagnostico |>
  select(serie, etiqueta_cat, pendiente, p_pendiente, potencia_130, lectura) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

cat("\nESTA TABLA Y LA DE PRECISIÓN VAN AL CAPÍTULO DE RESULTADOS.\n")
cat("Distinguen los tres casos que no se pueden confundir:\n")
cat("  - se detectó sesgo;\n")
cat("  - no se detectó y el diseño podía detectarlo (evidencia de ausencia);\n")
cat("  - no se detectó y el diseño no podía detectarlo (nada que concluir).\n")

# --- Redacción sugerida ---
cat("\nREDACCIÓN SUGERIDA\n")

n_informativo    <- sum(grepl("INFORMATIVO",    diagnostico$lectura))
n_no_concluyente <- sum(grepl("NO CONCLUYENTE", diagnostico$lectura))
n_detectado      <- sum(diagnostico$rechaza)

b_medio_HA       <- mean(flb_regresiones$pendiente[flb_regresiones$categoria != "D"])
agregado_rechaza <- any(flb_agregado$p_valor < 0.05)

# ZONA DE BORDE. Dicotomizar en 0,05 es exactamente lo que el resto del script
# desaconseja. Si el test agregado no cruza el umbral pero queda muy cerca y
# con la pendiente del lado del FLB, decir que el rechazo individual "no
# sobrevive a la agregación" es falso: p = 0,045 y p = 0,051 son el mismo
# resultado. Esa situación se trata aparte y ANTES que las demás.
p_agr_min      <- min(flb_agregado$p_valor)
agregado_borde <- !agregado_rechaza && p_agr_min < 0.10 &&
                  all(flb_agregado$pendiente > 1)

if (agregado_borde) {
  # Todo el texto de esta rama se construye con los valores estimados: no debe
  # quedar ningún número escrito a mano, porque una reejecución con datos
  # distintos lo dejaría desactualizado sin avisar.
  i_borde <- which.min(flb_agregado$p_valor)
  cat(sprintf("El test agregado queda en la zona de borde: p = %.4f, con b = %.4f\n",
              p_agr_min, flb_agregado$pendiente[i_borde]))
  cat(sprintf("e IC 95%% [%.3f ; %.3f], es decir del lado del FLB.\n",
              flb_agregado$ic95_inf[i_borde], flb_agregado$ic95_sup[i_borde]))
  cat("\nNO redactar esto como 'no se detectó sesgo' ni como 'el rechazo\n")
  cat("individual no sobrevive a la agregación': con p individual mínimo de\n")
  cat(sprintf("%.4f y p agregado de %.4f se trata del MISMO resultado, y la\n",
              min(flb_regresiones$p_pendiente), p_agr_min))
  cat("diferencia es de qué lado del umbral convencional cae, no de qué dice\n")
  cat("el dato.\n")
  cat("\nREDACCIÓN CORRECTA: encabezar con la ESTIMACIÓN y su INTERVALO, no\n")
  cat(sprintf("con la etiqueta de significación. 'La pendiente estimada es b = %.2f,\n",
              flb_agregado$pendiente[i_borde]))
  cat(sprintf("con un intervalo del 95%% de [%.3f ; %.3f], que apenas contiene el\n",
              flb_agregado$ic95_inf[i_borde], flb_agregado$ic95_sup[i_borde]))
  cat("valor de eficiencia.\n")
  cat("La evidencia apunta de forma consistente a un sesgo de tipo\n")
  cat("favourite-longshot de magnitud moderada, sin que el diseño alcance a\n")
  cat("separarlo de la ausencia de sesgo con el nivel de confianza\n")
  cat("convencional.' Ver la tabla de precisión para qué valores de b sí\n")
  cat("quedan descartados.\n")

} else if (n_detectado == 0 && agregado_rechaza) {
  cat("Ninguna regresión individual rechaza, pero el test agregado sí. Es el\n")
  cat("caso esperable cuando el efecto es real y moderado: las regresiones\n")
  cat("por separado carecen de precisión y la agregación la recupera.\n")
  cat("Encabezar con el test agregado y presentar las seis como desagregación.\n")

} else if (n_detectado > 0 && !agregado_rechaza) {
  cat("Alguna regresión individual rechaza pero el agregado no. Con seis\n")
  cat("contrastes eso es compatible con el azar. NO encabezar con el rechazo\n")
  cat("aislado: reportar el agregado como resultado principal y el rechazo\n")
  cat("individual como observación desagregada, señalando explícitamente que\n")
  cat("no sobrevive a la agregación.\n")

} else if (n_detectado == 0 && !agregado_rechaza) {
  cat("Ningún contraste rechaza la eficiencia. La formulación correcta\n")
  cat("depende de la potencia, y NO es homogénea entre categorías:\n")
  if (n_informativo > 0) {
    cat(sprintf("  - en %d casos el diseño tenía potencia adecuada: ahí sí se puede\n",
                n_informativo))
    cat("    afirmar que no hay FLB de magnitud relevante;\n")
  }
  if (n_no_concluyente > 0) {
    cat(sprintf("  - en %d casos no la tenía: ahí no hay nada que concluir y debe\n",
                n_no_concluyente))
    cat("    declararse como tal.\n")
  }
  cat(sprintf("\nAdemás: la pendiente media en H y A es %.3f. Si los intervalos\n", b_medio_HA))
  cat("incluyen valores sustancialmente mayores que 1, la conclusión honesta\n")
  cat("no es 'no hay sesgo' sino 'el estudio descarta un sesgo grande pero\n")
  cat("no tiene precisión para descartar uno leve'. Ver la tabla de precisión.\n")

} else {
  cat("Hay rechazos tanto individuales como en el test agregado. Encabezar\n")
  cat("con el agregado, reportar magnitud e intervalo, y verificar que ambas\n")
  cat("series coincidan antes de afirmar robustez.\n")
}

cat("\nEsta sección también justifica cuantitativamente el uso de cuatro\n")
cat("temporadas: con 380 partidos la potencia sería una fracción de la\n")
cat("obtenida acá, lo que convierte esa decisión en un argumento medido y no\n")
cat("en una preferencia.\n")



# -----------------------------------------------------------------------------
# 8. DIAGRAMA DE FIABILIDAD POR DECILES (VISUALIZACIÓN)
# -----------------------------------------------------------------------------
# Los deciles se conservan como GRÁFICO, no como inferencia. Cada punto es un
# decil de probabilidad implícita; el eje Y es la frecuencia real observada
# en ese decil, con intervalo de confianza binomial de Wilson.
#
# El intervalo de Wilson se prefiere al normal (Wald) porque se comporta mejor
# con proporciones cercanas a 0 o 1 y con tamaños de bin moderados, que es
# exactamente la situación de los deciles extremos.
#
# LIMITACIÓN DEL INTERVALO (declararla en el texto). Cada decil mezcla las tres
# categorías, de modo que un mismo partido puede aportar dos observaciones al
# mismo decil (por ejemplo, empate cotizado a 0,28 y visitante a 0,30). El
# intervalo de Wilson supone ensayos de Bernoulli independientes, supuesto que
# acá no se cumple exactamente. La consecuencia práctica es que los intervalos
# son algo más angostos de lo que corresponde y el recuento de "deciles fuera
# del IC" resulta levemente anti-conservador. Por eso este bloque es
# ESTRICTAMENTE DESCRIPTIVO: la inferencia queda en la Sección 7, que sí
# agrupa por partido.

wilson_inf <- function(exitos, n, nivel = 0.95) {
  z <- qnorm(1 - (1 - nivel) / 2)
  p <- exitos / n
  centro <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  radio  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  pmax(0, centro - radio)
}

wilson_sup <- function(exitos, n, nivel = 0.95) {
  z <- qnorm(1 - (1 - nivel) / 2)
  p <- exitos / n
  centro <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  radio  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  pmin(1, centro + radio)
}

construir_bins <- function(df, columna_pi, etiqueta_serie, n_bins = 10) {
  df |>
    filter(!is.na(.data[[columna_pi]])) |>
    mutate(bin = ntile(.data[[columna_pi]], n_bins)) |>
    summarise(
      n           = n(),
      pi_med      = mean(.data[[columna_pi]]),
      ocurrencias = sum(ocurrio),
      frec_real   = mean(ocurrio),
      .by = bin
    ) |>
    mutate(
      inf   = wilson_inf(ocurrencias, n),
      sup   = wilson_sup(ocurrencias, n),
      serie = etiqueta_serie,
      sesgo = pi_med - frec_real
    ) |>
    arrange(bin)
}

bins_mercado <- bind_rows(
  construir_bins(mercado_largo, "pi_precierre", "Pre-cierre"),
  construir_bins(mercado_largo, "pi_cierre",    "Cierre")
)

cat("\n=== DIAGRAMA DE FIABILIDAD POR DECILES (4 temporadas) ===\n")
bins_mercado |>
  select(serie, bin, n, pi_med, frec_real, inf, sup, sesgo) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf)

# --- Gráfico 5: diagrama de fiabilidad con IC ---
g_fiabilidad <- bins_mercado |>
  ggplot(aes(x = pi_med, y = frec_real, color = serie)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbar(aes(ymin = inf, ymax = sup), width = 0.012, alpha = 0.7) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("Pre-cierre" = "#2c3e50",
                                "Cierre"     = "#16a085")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Diagrama de fiabilidad del mercado",
    subtitle = "Deciles de probabilidad implícita con IC de Wilson al 95% | 2021/22–2024/25\nDiagonal = calibración perfecta",
    x        = "Probabilidad implícita normalizada (promedio del decil)",
    y        = "Frecuencia real de ocurrencia",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g_fiabilidad)

# --- Correlación descriptiva (continuidad con la versión anterior) ---
# Se conserva como estadístico DESCRIPTIVO. La inferencia la hace la Sección 7.
correlaciones <- bins_mercado |>
  summarise(correlacion = cor(pi_med, sesgo), .by = serie) |>
  mutate(correlacion = round(correlacion, 4))

cat("\n--- Correlación descriptiva (pi media del decil vs. sesgo del decil) ---\n")
print(correlaciones)
cat("Nota: estadístico descriptivo, sin error estándar ni test. Solo 10 puntos.\n")
cat("      La inferencia sobre FLB está en la Sección 7.\n")

# 8.1 CONCLUSIÓN DEL DIAGRAMA DE FIABILIDAD
# -----------------------------------------------------------------------------

cat("\n", strrep("-", 78), "\n", sep = "")
cat("CONCLUSIÓN — SECCIÓN 8\n")
cat(strrep("-", 78), "\n", sep = "")

# Un decil está mal calibrado de forma detectable si su probabilidad
# implícita media cae FUERA del intervalo de Wilson de la frecuencia real.
calibracion_bins <- bins_mercado |>
  mutate(fuera_ic = pi_med < inf | pi_med > sup) |>
  summarise(
    deciles       = n(),
    fuera_ic      = sum(fuera_ic),
    sobrestima    = sum(fuera_ic & pi_med > sup),   # implícita > frecuencia real
    subestima     = sum(fuera_ic & pi_med < inf),
    .by = serie
  )

cat("\n--- Deciles con calibración significativamente desviada ---\n")
print(calibracion_bins, width = Inf)
cat("(sobrestima = la probabilidad implícita supera la frecuencia observada)\n")

for (i in seq_len(nrow(calibracion_bins))) {
  f <- calibracion_bins[i, ]
  cat(sprintf("\n%s: %d de %d deciles fuera del IC del 95%%.\n",
              f$serie, f$fuera_ic, f$deciles))
  if (f$fuera_ic == 0) {
    cat("  Ningún decil se aparta de forma detectable de la calibración perfecta.\n")
  } else if (f$sobrestima > 0 && f$subestima == 0) {
    cat("  Todas las desviaciones son sobrestimaciones del mercado.\n")
  } else if (f$subestima > 0 && f$sobrestima == 0) {
    cat("  Todas las desviaciones son subestimaciones del mercado.\n")
  } else {
    cat("  Desviaciones en ambos sentidos: no hay un patrón direccional único.\n")
  }
}

cat("\nCAUTELA: con 10 deciles y un IC del 95%, esperar cerca de un decil\n")
cat("fuera del intervalo por azar es lo normal aun con calibración perfecta.\n")
cat("Además, los IC de Wilson suponen observaciones independientes y cada\n")
cat("decil puede contener dos resultados del mismo partido, de modo que son\n")
cat("algo más angostos de lo que corresponde: el recuento es levemente\n")
cat("anti-conservador y NO debe usarse como contraste formal.\n")
cat("Lo relevante no es el recuento sino si las desviaciones siguen el patrón\n")
cat("del FLB: sobrestimación concentrada en los deciles bajos (longshots) y\n")
cat("subestimación en los altos (favoritos).\n")

# Patrón direccional: ¿el sesgo cae al subir la probabilidad implícita?
cat("\n--- Correlación descriptiva (estadístico, NO test) ---\n")
print(correlaciones)
for (i in seq_len(nrow(correlaciones))) {
  r <- correlaciones$correlacion[i]
  cat(sprintf("%s: r = %+.4f -> %s\n",
              correlaciones$serie[i], r,
              if (r < -0.5) "patrón compatible con FLB clásico"
              else if (r < 0) "patrón débilmente negativo, no concluyente"
              else "sin patrón negativo: no sugiere FLB"))
}
cat("\nRecordatorio: son 10 puntos sin error estándar. La inferencia está en\n")
cat("la Sección 7. Reportar esto como descriptivo y nada más.\n")



# -----------------------------------------------------------------------------
# 9. RETORNOS POR BANDA DE PRECIO (Cain, Law y Peel, 2000)
# -----------------------------------------------------------------------------
# 9.1 TABLA DE RETORNOS POR BANDA
# -----------------------------------------------------------------------------
# Replicación del test que Cain, Law y Peel (2000) aplican efectivamente a los
# resultados simples (su Tabla 1): el retorno medio a una apuesta unitaria,
# desagregado por banda de precio.
#
# Notación del paper: el "precio" es i = 1/(1 + a), con a la cuota en formato
# fraccionario (odds against). En cuota decimal o, eso es exactamente
# i = 1/o, es decir, la probabilidad implícita CRUDA (sin normalizar).
# Acá se usa esa definición para que las bandas sean directamente comparables
# con las suyas.
#
# Retorno de una apuesta unitaria a cuota decimal o:
#     +(o − 1)  si el resultado ocurre
#     −1        si no ocurre
#
# El FLB se manifiesta como retornos MENOS negativos en las bandas de precio
# alto (favoritos) que en las de precio bajo (longshots).
#
# Este análisis es complementario y de naturaleza económica: mide el sesgo en
# unidades monetarias, no probabilísticas. La inferencia formal sigue siendo
# la de la Sección 7.

# --- Estimador de varianza agrupado para lm (se usa en toda la Sección 9) ---
# Difiere del de la Sección 7.2 porque vcov() de un lm incluye el factor s^2,
# que no corresponde acá: el "pan" del sándwich debe ser (X'X)^-1 puro.
#
# Se define acá arriba, y no dentro del test de tendencia, porque TODOS los
# errores estándar de esta sección se agrupan por partido: la tabla por banda,
# el contraste entre extremos y la tendencia. Las tres apuestas de un mismo partido no son
# independientes, y en las bandas extremas esa dependencia es fuerte.
vcov_agrupado_lm <- function(modelo, grupo) {
  X <- model.matrix(modelo)
  e <- residuals(modelo)
  U <- rowsum(X * e, factor(grupo))
  pan   <- solve(crossprod(X))
  carne <- crossprod(U)
  G <- nlevels(factor(grupo)); n <- nrow(X); k <- ncol(X)
  correccion <- (G / (G - 1)) * ((n - 1) / (n - k))
  pan %*% carne %*% pan * correccion
}

retornos_por_banda <- function(df, columna_cuota, etiqueta_serie) {

  d <- df |>
    filter(!is.na(.data[[columna_cuota]])) |>
    mutate(
      id_partido = paste(date, home_team, away_team),
      precio     = 1 / .data[[columna_cuota]],        # i = 1/o (prob. cruda)
      retorno    = if_else(ocurrio == 1, .data[[columna_cuota]] - 1, -1),
      banda      = cut(precio,
                       breaks = c(0, 0.2, 0.4, 0.6, 1),
                       labels = c("0 < i ≤ 0,2", "0,2 < i ≤ 0,4",
                                  "0,4 < i ≤ 0,6", "0,6 < i ≤ 1"),
                       include.lowest = TRUE)
    )

  niveles <- levels(d$banda)

  # El error estándar de la media de cada banda se obtiene como el de la
  # constante de una regresión sin regresores, con varianza agrupada por
  # partido. Se conserva también el EE ingenuo para dejar visible el efecto
  # de la corrección.
  map_dfr(niveles, function(bb) {
    sub <- d |> filter(banda == bb)
    m   <- lm(retorno ~ 1, data = sub)
    tibble(
      banda         = factor(bb, levels = niveles),
      n             = nrow(sub),
      partidos      = n_distinct(sub$id_partido),
      precio_medio  = mean(sub$precio),
      retorno_medio = mean(sub$retorno),
      ee_ingenuo    = sd(sub$retorno) / sqrt(nrow(sub)),
      ee_retorno    = sqrt(vcov_agrupado_lm(m, sub$id_partido)[1, 1])
    )
  }) |>
    mutate(
      serie    = etiqueta_serie,
      ic95_inf = retorno_medio - qnorm(0.975) * ee_retorno,
      ic95_sup = retorno_medio + qnorm(0.975) * ee_retorno
    ) |>
    arrange(banda)
}

tabla_cain <- bind_rows(
  retornos_por_banda(mercado_largo, "cuota_precierre", "Pre-cierre"),
  retornos_por_banda(mercado_largo, "cuota_cierre",    "Cierre")
)

cat("\n=== RETORNOS A UNA APUESTA UNITARIA POR BANDA DE PRECIO ===\n")
cat("(replicación de Cain, Law y Peel, 2000, Tabla 1)\n")
cat("Precio i = 1/cuota decimal = probabilidad implícita cruda\n\n")
cat("Errores estándar agrupados por partido (ee_retorno); ee_ingenuo es el\n")
cat("que saldría de tratar las filas como independientes.\n\n")
tabla_cain |>
  select(serie, banda, n, partidos, precio_medio, retorno_medio,
         ee_ingenuo, ee_retorno, ic95_inf, ic95_sup) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

cat("\nInterpretación: si hay FLB, el retorno medio MEJORA (es menos negativo)\n")
cat("                a medida que sube el precio, es decir, hacia los favoritos.\n")

cat("\n--- Efecto de agrupar por partido dentro de cada banda ---\n")
tabla_cain |>
  transmute(serie, banda, cociente = round(ee_retorno / ee_ingenuo, 3)) |>
  print(n = Inf, width = Inf)
cat("Dentro de una misma banda los retornos de un mismo partido están\n")
cat("correlacionados NEGATIVAMENTE (a lo sumo uno de los tres resultados\n")
cat("ocurre), de modo que el cociente suele ser MENOR que 1: acá el EE\n")
cat("ingenuo era conservador. En el contraste entre bandas extremas de la\n")
cat("Sección 9.4 el efecto va en sentido contrario y sí importa.\n")



# -----------------------------------------------------------------------------
# 9.2 REFERENCIA: RETORNO ESPERADO BAJO CALIBRACIÓN PERFECTA
# -----------------------------------------------------------------------------
# La referencia natural del gráfico parecería ser cero, pero con overround
# ninguna apuesta puede rendir cero: la referencia correcta no es esa.
#
# Bajo calibración perfecta y margen uniforme entre los tres resultados,
# CUALQUIER apuesta unitaria tiene el mismo retorno esperado, sin importar
# la cuota:
#
#   E[R] = p*o - 1 = pi_norm / (pi_norm*(1+OR)) - 1 = 1/(1+OR) - 1
#
# Con un overround del 5,4%, eso da aproximadamente -5,1%. Las bandas por
# encima de esa línea reciben MENOS margen que el promedio; las que quedan
# por debajo, MÁS. Esa desviación según el nivel de precio es exactamente el
# favourite-longshot bias medido en unidades monetarias.

or_medio_serie <- c(
  "Pre-cierre" = mean(datos_mercado$overround,   na.rm = TRUE),
  "Cierre"     = mean(datos_mercado$overround_c, na.rm = TRUE)
)

referencia_justa <- tibble(
  serie         = names(or_medio_serie),
  overround     = as.numeric(or_medio_serie),
  retorno_justo = 1 / (1 + as.numeric(or_medio_serie)) - 1
)

cat("\n--- Retorno esperado bajo calibración perfecta y margen uniforme ---\n")
referencia_justa |>
  mutate(overround     = round(overround * 100, 3),
         retorno_justo = round(retorno_justo, 4)) |>
  print(width = Inf)

tabla_cain_ref <- tabla_cain |>
  left_join(referencia_justa |> select(serie, retorno_justo), by = "serie") |>
  mutate(
    exceso   = retorno_medio - retorno_justo,
    z_exceso = exceso / ee_retorno,
    p_exceso = 2 * pnorm(-abs(z_exceso)),
    lectura  = case_when(
      p_exceso >= 0.05 ~ "Compatible con margen uniforme",
      exceso   <  0    ~ "Margen MAYOR al promedio (peor que justo)",
      TRUE             ~ "Margen MENOR al promedio (mejor que justo)"
    )
  )

cat("\n--- Exceso de retorno sobre la referencia justa, por banda ---\n")
tabla_cain_ref |>
  select(serie, banda, n, retorno_medio, retorno_justo, exceso,
         z_exceso, p_exceso, lectura) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

cat("\nLECTURA. Exceso negativo en las bandas de precio bajo (longshots) y\n")
cat("positivo en las de precio alto (favoritos) es la firma del FLB medida\n")
cat("en pesos: la casa carga más margen sobre los longshots.\n")


# -----------------------------------------------------------------------------
# 9.3 GRÁFICO: RETORNOS POR BANDA CON LA REFERENCIA CORRECTA
# -----------------------------------------------------------------------------

g_cain <- tabla_cain_ref |>
  ggplot(aes(x = banda, y = retorno_medio, color = serie, group = serie)) +
  geom_hline(aes(yintercept = retorno_justo, color = serie),
             linetype = "dashed", linewidth = 0.7, alpha = 0.8) +
  geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup),
                  position = position_dodge(width = 0.3), size = 0.6) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("Pre-cierre" = "#2c3e50",
                                "Cierre"     = "#16a085")) +
  labs(
    title    = "Retorno a una apuesta unitaria por banda de precio",
    subtitle = "Réplica de Cain, Law y Peel (2000, Tabla 1) | Premier League 2021/22–2024/25\nIC del 95% | Línea discontinua: retorno esperado bajo calibración perfecta y margen uniforme",
    x        = "Banda de precio i = 1/cuota",
    y        = "Retorno medio por unidad apostada",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g_cain)



# -----------------------------------------------------------------------------
# 9.4 CONCLUSIÓN DEL TEST DE RETORNOS POR BANDA
# -----------------------------------------------------------------------------

cat("\n", strrep("-", 78), "\n", sep = "")
cat("CONCLUSIÓN — SECCIÓN 9\n")
cat(strrep("-", 78), "\n", sep = "")

# --- Contraste entre la banda de precio más baja y la más alta ---
# CORRECCIÓN (versión 4). La versión anterior trataba las dos medias como
# INDEPENDIENTES. No lo son, y no por un tecnicismo: las dos bandas extremas
# están pobladas por los MISMOS partidos. Un partido con un favorito a 1,40
# (i ≈ 0,71, banda alta) tiene por construcción un visitante alrededor de 7,00
# (i ≈ 0,14, banda baja). Y dentro del partido los retornos están fuertemente
# correlacionados en sentido NEGATIVO: si gana el favorito, R_alta = +0,40 y
# R_baja = −1; si gana el visitante, R_alta = −1 y R_baja = +6.
#
# Como Var(R_alta − R_baja) = Var_a + Var_b − 2·Cov y acá Cov < 0, el error
# estándar ingenuo SUBESTIMA la incertidumbre y el contraste resultaba
# anti-conservador. La corrección: estimar la diferencia como el coeficiente
# de una variable indicadora de banda alta sobre las observaciones
# individuales de las dos bandas extremas, con varianza agrupada por partido.
# Es el mismo criterio de la Sección 7.2 y del test de tendencia.
#
# Bajo FLB se espera R_alta > R_baja, es decir diferencia > 0.

contraste_extremos_fn <- function(df, columna_cuota, etiqueta_serie) {

  d <- df |>
    filter(!is.na(.data[[columna_cuota]])) |>
    mutate(
      id_partido = paste(date, home_team, away_team),
      precio     = 1 / .data[[columna_cuota]],
      retorno    = if_else(ocurrio == 1, .data[[columna_cuota]] - 1, -1)
    ) |>
    filter(precio <= 0.2 | precio > 0.6) |>
    mutate(alta = as.integer(precio > 0.6))

  m      <- lm(retorno ~ alta, data = d)
  V_cr   <- vcov_agrupado_lm(m, d$id_partido)
  dif    <- unname(coef(m)["alta"])
  se_cr  <- sqrt(V_cr["alta", "alta"])
  # EE "ingenuo": el que daba la versión anterior, tratando las dos medias
  # como muestras independientes.
  se_ing <- sqrt(var(d$retorno[d$alta == 1]) / sum(d$alta == 1) +
                 var(d$retorno[d$alta == 0]) / sum(d$alta == 0))

  # Cuántos partidos aportan observaciones a las DOS bandas a la vez: es la
  # medida directa de por qué el supuesto de independencia no se sostiene.
  en_ambas <- tapply(d$alta, d$id_partido, function(x) length(unique(x)) == 2)

  z <- dif / se_cr

  tibble(
    serie          = etiqueta_serie,
    banda_baja     = "0 < i <= 0,2",
    banda_alta     = "0,6 < i <= 1",
    n_baja         = sum(d$alta == 0),
    n_alta         = sum(d$alta == 1),
    partidos       = n_distinct(d$id_partido),
    partidos_ambas = sum(en_ambas),
    r_baja         = mean(d$retorno[d$alta == 0]),
    r_alta         = mean(d$retorno[d$alta == 1]),
    diferencia     = dif,
    ee_ingenuo     = se_ing,
    ee_dif         = se_cr,
    z              = z,
    p_valor        = 2 * pnorm(-abs(z))
  )
}

contraste_extremos <- bind_rows(
  contraste_extremos_fn(mercado_largo, "cuota_precierre", "Pre-cierre"),
  contraste_extremos_fn(mercado_largo, "cuota_cierre",    "Cierre")
)

cat("\n--- Contraste entre banda de precio más alta y más baja ---\n")
cat("Diferencia estimada sobre observaciones individuales, EE agrupado por partido\n\n")
contraste_extremos |>
  select(serie, n_baja, n_alta, partidos, partidos_ambas,
         r_baja, r_alta, diferencia, ee_ingenuo, ee_dif, z, p_valor) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(width = Inf)

cat("\npartidos_ambas es la cantidad de partidos que aportan observaciones a\n")
cat("las DOS bandas simultáneamente: es la razón por la que las dos medias no\n")
cat("pueden tratarse como muestras independientes.\n")

for (i in seq_len(nrow(contraste_extremos))) {
  f <- contraste_extremos[i, ]
  cat(sprintf("\n%s:\n", f$serie))
  cat(sprintf("  Retorno en longshots: %+.4f | en favoritos: %+.4f | dif: %+.4f\n",
              f$r_baja, f$r_alta, f$diferencia))
  cat(sprintf("  EE ingenuo: %.4f | EE agrupado: %.4f (cociente %.3f)\n",
              f$ee_ingenuo, f$ee_dif, f$ee_dif / f$ee_ingenuo))
  if (f$ee_dif > f$ee_ingenuo) {
    cat("  El agrupamiento AMPLÍA el error estándar, como se esperaba: los\n")
    cat("  retornos de un mismo partido en las dos bandas están correlacionados\n")
    cat("  negativamente. El contraste sin agrupar era anti-conservador.\n")
  }
  if (f$p_valor < 0.05 && f$diferencia > 0) {
    cat("  Los favoritos rinden significativamente mejor que los longshots:\n")
    cat("  evidencia de FLB en términos económicos.\n")
  } else if (f$p_valor < 0.05 && f$diferencia < 0) {
    cat("  Los longshots rinden significativamente mejor: patrón inverso al FLB.\n")
  } else {
    cat(sprintf("  Diferencia en la dirección del FLB pero no significativa (p = %.4f).\n",
                f$p_valor))
  }
}

cat("\nNOTA. Este contraste usa solo las dos bandas extremas y descarta las\n")
cat("intermedias. El test de tendencia que sigue usa las 4.560 observaciones\n")
cat("y es el instrumento económico de referencia; este queda como el número\n")
cat("directamente comparable con Cain, Law y Peel (2000).\n")

# --- Monotonía a lo largo de las cuatro bandas (DESCRIPTIVO) ---
# Se conserva como descriptivo, pero NO como evidencia formal: exigir
# monotonía estricta sobre cuatro medias es un criterio frágil. Dos bandas
# separadas por menos de un error estándar pueden invertirse por azar y
# romper la condición sin que el patrón subyacente cambie. El contraste
# formal de la tendencia está en el bloque siguiente.
cat("\n--- Monotonía del retorno a lo largo de las bandas (descriptivo) ---\n")
monotonia <- tabla_cain |>
  arrange(serie, banda) |>
  summarise(
    creciente   = all(diff(retorno_medio) > 0),
    decreciente = all(diff(retorno_medio) < 0),
    salto_min   = min(abs(diff(retorno_medio))),
    .by = serie
  )
print(monotonia, width = Inf)
cat("El FLB predice retorno creciente con el precio (hacia los favoritos).\n")
cat("ATENCIÓN a la columna salto_min: si el menor salto entre bandas es\n")
cat("mucho más chico que los errores estándar de la tabla anterior, esas dos\n")
cat("bandas son indistinguibles y la monotonía estricta no significa nada.\n")
cat("Por eso el contraste formal lo hace el test de tendencia que sigue.\n")


# -----------------------------------------------------------------------------
# TEST DE TENDENCIA DEL RETORNO SEGÚN EL PRECIO
# -----------------------------------------------------------------------------
# Reemplaza al chequeo de monotonía como evidencia formal, y también mejora
# al contraste entre bandas extremas: en vez de comparar dos medias, usa las
# 4.560 observaciones individuales.
#
#     retorno_i = alfa + beta * precio_i + e_i
#
# Bajo calibración perfecta y margen uniforme, el retorno esperado es
# constante e igual a 1/(1+OR) − 1 para cualquier precio, de modo que
# beta = 0. El FLB implica beta > 0: el retorno mejora hacia los favoritos.
#
# Los errores estándar se agrupan por partido, por la misma razón que en la
# Sección 7.2: las tres apuestas de un mismo partido no son independientes.

# El estimador de varianza agrupado para lm (vcov_agrupado_lm) está definido
# al comienzo de la Sección 9: desde la versión 4 lo usan también la tabla por
# banda y el contraste entre extremos.

tendencia_retornos <- function(df, columna_cuota, etiqueta_serie) {

  d <- df |>
    mutate(id_partido = paste(date, home_team, away_team)) |>
    filter(!is.na(.data[[columna_cuota]])) |>
    transmute(
      id_partido,
      precio  = 1 / .data[[columna_cuota]],
      retorno = if_else(ocurrio == 1, .data[[columna_cuota]] - 1, -1)
    )

  m    <- lm(retorno ~ precio, data = d)
  V_cr <- vcov_agrupado_lm(m, d$id_partido)

  beta <- unname(coef(m)["precio"])
  se   <- sqrt(V_cr["precio", "precio"])
  z    <- beta / se

  tibble(
    serie         = etiqueta_serie,
    observaciones = nrow(d),
    partidos      = n_distinct(d$id_partido),
    beta          = beta,
    se_agrupado   = se,
    ic95_inf      = beta - qnorm(0.975) * se,
    ic95_sup      = beta + qnorm(0.975) * se,
    z             = z,
    p_valor       = 2 * pnorm(-abs(z))
  )
}

tendencia_cain <- bind_rows(
  tendencia_retornos(mercado_largo, "cuota_precierre", "Pre-cierre"),
  tendencia_retornos(mercado_largo, "cuota_cierre",    "Cierre")
)

cat("\n--- Test de tendencia: retorno = alfa + beta * precio ---\n")
cat("H0: beta = 0 (margen uniforme) | FLB implica beta > 0\n")
cat("Errores estándar agrupados por partido | 4.560 observaciones\n\n")

tendencia_cain |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(width = Inf)

for (i in seq_len(nrow(tendencia_cain))) {
  f <- tendencia_cain[i, ]
  cat(sprintf("\n%s: beta = %+.4f | IC 95%% [%+.4f ; %+.4f] | p = %.4f\n",
              f$serie, f$beta, f$ic95_inf, f$ic95_sup, f$p_valor))
  if (f$p_valor < 0.05 && f$beta > 0) {
    cat("  El retorno crece significativamente con el precio: evidencia de FLB\n")
    cat("  en términos económicos, contrastada sobre las observaciones\n")
    cat("  individuales y no sobre cuatro medias.\n")
  } else if (f$p_valor < 0.05) {
    cat("  El retorno DECRECE con el precio: patrón inverso al FLB.\n")
  } else {
    cat("  Tendencia en la dirección del FLB pero no significativa.\n")
  }
  cat(sprintf("  Magnitud: pasar de un precio de 0,15 (longshot) a 0,70\n"))
  cat(sprintf("  (favorito) mejora el retorno esperado en %+.4f unidades por\n",
              f$beta * 0.55))
  cat("  unidad apostada.\n")
}

# --- Potencia de este test frente a la regresión ---
cat("\n--- Por qué la regresión es el instrumento primario ---\n")
cat("Los retornos de apuestas tienen varianza enorme: un longshot paga +6 o\n")
cat("-1, sin valores intermedios. Por eso los errores estándar de esta tabla\n")
cat("son grandes y el test económico es MENOS eficiente que la regresión de\n")
cat("calibración, que aprovecha toda la información contenida en la\n")
cat("probabilidad y no solo el pago realizado.\n")
cat("Conviene decirlo explícitamente en el texto: esta tabla es el resultado\n")
cat("comparable con Cain, Law y Peel (2000), y la regresión es el que sostiene\n")
cat("la inferencia.\n")

# --- Coherencia con la Sección 7 ---
# Se comparan DIRECCIÓN y MAGNITUD, no banderas de significación. Comparar
# solo si cada p-valor cruza 0,05 convierte una diferencia entre p = 0,045 y
# p = 0,077 en una contradicción, cuando en realidad ambos apuntan al mismo
# lado con magnitudes similares.
cat("\n--- Coherencia con el test de la Sección 7 ---\n")

dir_regresion <- mean(flb_regresiones$pendiente[flb_regresiones$categoria != "D"]) - 1
dir_economico <- mean(contraste_extremos$diferencia)
dir_tendencia <- mean(tendencia_cain$beta)

p_min_regresion <- min(flb_regresiones$p_pendiente)
p_min_agregado  <- min(flb_agregado$p_valor)
# El lado económico tiene ahora dos contrastes: el de bandas extremas (dos
# medias) y el de tendencia (4.560 observaciones). Se toma el mínimo, pero se
# reportan por separado porque el de tendencia es el más potente de los dos.
p_min_extremos  <- min(contraste_extremos$p_valor)
p_min_tendencia <- min(tendencia_cain$p_valor)
p_min_economico <- min(p_min_extremos, p_min_tendencia)

cat(sprintf("Regresión (H y A): pendiente media = %.4f, es decir %+.4f respecto de 1\n",
            dir_regresion + 1, dir_regresion))
cat(sprintf("                   p mínimo individual = %.4f | p del test agregado = %.4f\n",
            p_min_regresion, p_min_agregado))
cat(sprintf("Económico:         diferencia media favoritos - longshots = %+.4f\n",
            dir_economico))
cat(sprintf("                   p bandas extremas = %.4f | p tendencia = %.4f\n",
            p_min_extremos, p_min_tendencia))
cat(sprintf("Tendencia:         beta medio = %+.4f (positivo = consistente con FLB)\n",
            mean(tendencia_cain$beta)))

# Se exige coincidencia de signo en los TRES enfoques, no en dos
misma_direccion <- (dir_regresion > 0) == (dir_economico > 0) &&
                   (dir_regresion > 0) == (dir_tendencia > 0)

cat("\n")
if (misma_direccion) {
  cat("Los dos enfoques apuntan en la MISMA dirección.\n")
  if (p_min_agregado < 0.05 && p_min_economico < 0.05) {
    cat("Ambos alcanzan significación: la evidencia de FLB es consistente en\n")
    cat("escala probabilística y en retorno económico.\n")
  } else if (p_min_agregado < 0.05 || p_min_economico < 0.05) {
    cat("Uno alcanza significación y el otro queda al margen, pero el signo\n")
    cat("coincide. NO es una contradicción: es lo esperable cuando el efecto\n")
    cat("es real pero pequeño respecto de la precisión disponible, y el test\n")
    cat("económico tiene menos potencia por la varianza de los pagos.\n")
    cat("Redactar como evidencia CONVERGENTE PERO NO CONCLUYENTE.\n")
  } else {
    cat("Ninguno alcanza significación, pero ambos apuntan en la misma\n")
    cat("dirección. Reportar la consistencia de signo y la falta de precisión,\n")
    cat("sin presentarlo como confirmación del sesgo ni como su ausencia.\n")
  }
} else {
  cat("ATENCIÓN: los dos enfoques apuntan en direcciones OPUESTAS.\n")
  cat("Causa más probable: el test de retornos agrega las tres categorías y\n")
  cat("usa precios crudos sin normalizar, de modo que mezcla el efecto del\n")
  cat("overround con el del sesgo. La regresión, que trabaja con\n")
  cat("probabilidades normalizadas y por categoría, es la referencia\n")
  cat("inferencial. Reportar la discrepancia y explicarla.\n")
}

cat("\nRecordatorio: todos los retornos son negativos por construcción, dado\n")
cat("que el overround garantiza valor esperado negativo. Lo que se contrasta\n")
cat("es el PATRÓN entre bandas respecto de la referencia justa, no el nivel.\n")



# -----------------------------------------------------------------------------
# 10. OBJETO DE VALIDACIÓN 2024/25 PARA LOS SCRIPTS 08 Y 09
# -----------------------------------------------------------------------------
# A partir de acá se cambia de ventana temporal: la comparación modelo vs.
# mercado involucra al modelo y por lo tanto DEBE ser fuera de muestra.
# Solo se usan los 380 partidos de 2024/25.
#
# Los objetos datos, datos_largo y datos_bins son los que leen los Scripts 08
# y 09, e incluyen las columnas de ambas series de cuotas.

datos <- datos_validacion

datos_largo <- datos |>
  select(home_team, away_team, date, ftr,
         pi_H,   pi_D,   pi_A,
         piC_H,  piC_D,  piC_A,
         prob_H, prob_D, prob_A) |>
  pivot_longer(
    cols      = c(pi_H, pi_D, pi_A, piC_H, piC_D, piC_A, prob_H, prob_D, prob_A),
    names_to  = c(".value", "resultado"),
    names_sep = "_"
  ) |>
  rename(pi_mercado        = pi,
         pi_mercado_cierre = piC,
         prob_modelo       = prob) |>
  mutate(ocurrio = as.integer(ftr == resultado))

cat("\n=== OBJETO DE VALIDACIÓN (Scripts 08 y 09) ===\n")
cat("Partidos:", nrow(datos), "| Filas en formato largo:", nrow(datos_largo), "\n")

n_sin_cierre_val <- sum(is.na(datos$piC_H))
cat("Partidos de 24/25 sin cuota de cierre:", n_sin_cierre_val, "\n")
if (n_sin_cierre_val > 0) {
  cat("  -> El Script 08 calculará las métricas de la serie de cierre sobre\n")
  cat("     los", nrow(datos) - n_sin_cierre_val, "partidos completos, y lo declarará explícitamente.\n")
}

# datos_bins: deciles sobre la validación con la serie pre-cierre.
# Construcción IDÉNTICA a la versión anterior (necesaria para el chequeo de
# verificación 4 y para que el Script 09 no cambie).
n_bins <- 10

datos_bins <- datos_largo |>
  mutate(bin = ntile(pi_mercado, n_bins)) |>
  group_by(bin) |>
  summarise(
    n               = n(),
    pi_mercado_med  = mean(pi_mercado),
    prob_modelo_med = mean(prob_modelo),
    frec_real       = mean(ocurrio),
    .groups = "drop"
  ) |>
  mutate(
    sesgo_mercado = pi_mercado_med - frec_real,
    sesgo_modelo  = prob_modelo_med - frec_real
  )

cat("\n--- Deciles sobre la validación 24/25 (serie pre-cierre) ---\n")
datos_bins |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = 10)

# -----------------------------------------------------------------------------
# 10.1 ADVERTENCIA PARA LOS SCRIPTS 08 Y 09
# -----------------------------------------------------------------------------
# El patrón de calibración por deciles de la temporada de validación NO es el
# mismo que el de las cuatro temporadas. La correlación entre la probabilidad
# implícita media del decil y su sesgo cambia de SIGNO. No es un error de
# cálculo: son 380 partidos resumidos en 10 puntos, muestra en la que el
# patrón se pierde en el ruido. Pero tiene una consecuencia directa sobre cómo
# se redactan los Scripts 08 y 09, así que conviene dejarlo dicho acá.

cor_val_mercado <- cor(datos_bins$pi_mercado_med, datos_bins$sesgo_mercado)
cor_4temp_pre   <- correlaciones$correlacion[correlaciones$serie == "Pre-cierre"]

cat("\n--- ADVERTENCIA PARA LOS SCRIPTS 08 Y 09 ---\n")
cat(sprintf("Correlación pi-sesgo por decil, 4 temporadas (pre-cierre): %+.4f\n",
            cor_4temp_pre))
cat(sprintf("Correlación pi-sesgo por decil, solo 2024/25:              %+.4f\n",
            cor_val_mercado))

if (sign(cor_4temp_pre) != sign(cor_val_mercado)) {
  cat("\nEL SIGNO SE INVIERTE. Sobre las cuatro temporadas el patrón de deciles\n")
  cat("es compatible con el FLB clásico; sobre los 380 partidos de la temporada\n")
  cat("de validación, no lo es.\n")
  cat("\nCONSECUENCIA: la comparación modelo vs. mercado de los Scripts 08 y 09\n")
  cat("se hace sobre una temporada que NO exhibe el sesgo documentado en la\n")
  cat("Sección 7. Hay que declararlo explícitamente ahí. El mercado de 2024/25\n")
  cat("no es una versión en miniatura del mercado del período, y una eventual\n")
  cat("ventaja o desventaja del modelo frente al mercado en esa temporada no\n")
  cat("se puede extrapolar hacia atrás. La dirección correcta de lectura es la\n")
  cat("inversa: el estudio del mercado usa cuatro temporadas justamente porque\n")
  cat("una sola no alcanza para caracterizarlo.\n")
} else {
  cat("\nEl signo coincide con el del período completo: la temporada de\n")
  cat("validación reproduce cualitativamente el patrón general.\n")
}


# -----------------------------------------------------------------------------
# 11. CHEQUEOS DE VERIFICACIÓN OBLIGATORIOS
# -----------------------------------------------------------------------------
# Si alguno de estos no coincide, algo se rompió en el Script 01.
# (Los chequeos 1 y 2 —cantidad de filas y parámetros de Maher— están en el
#  Script 01 y en la reejecución de los Scripts 02 a 06 respectivamente.)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("CHEQUEOS DE VERIFICACIÓN\n")
cat(strrep("=", 78), "\n", sep = "")

# --- Chequeo 3: overround pre-cierre sobre 24/25 ≈ 5,5% ---
overround_2425 <- mean(datos$overround) * 100
ok_3 <- abs(overround_2425 - 5.5) < 0.5
cat(sprintf("3. Overround pre-cierre 24/25: %.4f%% (esperado ≈ 5,5%%) -> %s\n",
            overround_2425, if (ok_3) "OK" else "REVISAR"))

# --- Chequeo 4: correlación de sesgo sobre 24/25 con pre-cierre ≈ +0,11 ---
cor_2425 <- cor(datos_bins$pi_mercado_med, datos_bins$sesgo_mercado)
ok_4 <- abs(cor_2425 - 0.11) < 0.05
cat(sprintf("4. Correlación de sesgo 24/25:  %+.4f (esperado ≈ +0,11) -> %s\n",
            cor_2425, if (ok_4) "OK" else "REVISAR"))

if (!ok_3 || !ok_4) {
  warning("Al menos un chequeo de verificación falló. ",
          "Revisar el Script 01 antes de dar por buenos los cambios.")
} else {
  cat("\nAmbos chequeos del Script 07 pasan: la serie pre-cierre reproduce\n")
  cat("los resultados de la versión anterior.\n")
}


# -----------------------------------------------------------------------------
# 12. GUARDAR
# -----------------------------------------------------------------------------

saveRDS(
  list(
    # --- Validación 24/25: insumos de los Scripts 08 y 09 ---
    datos               = datos,
    datos_largo         = datos_largo,
    datos_bins          = datos_bins,

    # --- Análisis de mercado: 4 temporadas ---
    datos_mercado       = datos_mercado,
    mercado_largo       = mercado_largo,
    bins_mercado        = bins_mercado,

    # --- Resultados del análisis de mercado ---
    overround_global    = overround_global,
    overround_temporada = overround_temporada,
    pareado_temporada   = pareado_temporada,
    test_t_pareado      = t_pareado,
    test_w_pareado      = w_pareado,
    flb_regresiones     = flb_regresiones,
    flb_agregado        = flb_agregado,
    precision_flb       = precision,
    diagnostico_flb     = diagnostico,
    tabla_potencia      = tabla_potencia,
    correlaciones       = correlaciones,
    tabla_cain          = tabla_cain,
    tabla_cain_ref      = tabla_cain_ref,
    referencia_justa    = referencia_justa,
    contraste_extremos  = contraste_extremos,
    tendencia_cain      = tendencia_cain,

    # --- Diagnóstico de calidad de datos ---
    overround_extremos  = extremos,
    robustez_pareado    = if (hay_extremos) comparacion_robustez else NULL,
    robustez_temporal   = robustez_temporal
  ),
  "resultados_mercado.rds"
)

cat("\nArchivo 'resultados_mercado.rds' guardado → usar en Scripts 08 y 09.\n")
