# =============================================================================
# TESINA - Agustín Valente (904158)
# Análisis actuarial de cuotas en apuestas deportivas
# =============================================================================
# SCRIPT 08: Calibración probabilística
# =============================================================================
# Depende de: Script 07 (resultados_mercado.rds)
#
# Objetivos:
#   1. Calcular Brier Score global — mercado vs. modelo (Brier, 1950)
#   2. Calcular log-loss global  — mercado vs. modelo (Hvattum & Arntzen, 2010)
#   3. Contrastar formalmente si las diferencias entre competidores son
#      distinguibles de cero, mediante tests pareados sobre los scores
#      partido a partido (Diebold & Mariano, 1995)
#   4. Descomponer el Brier Score en calibración, resolución e incertidumbre
#      (Murphy, 1973) para separar las fuentes de error
#   5. Comparar ambas métricas por tipo de resultado (H, D, A)
#   6. Comparar ambas métricas por decil de probabilidad
#   7. Producir tabla resumen y gráficos publicables
#
# -----------------------------------------------------------------------------
# DECISIONES DE DISEÑO
#
#   (a) BRIER SCORE GLOBAL. Se calcula DIRECTAMENTE desde datos (formato
#       ancho, 380 filas), implementando literalmente la ecuación (2) de
#       Brier (1950):
#
#         BS = (1/n) * Σᵢ Σⱼ (fᵢⱼ − oᵢⱼ)²
#
#       donde la suma interna recorre los r = 3 resultados de cada partido i.
#       Rango: [0, 2]. Promediar el error cuadrático fila a fila sobre las
#       1.140 filas de datos_largo daría la fórmula clásica DIVIDIDA por r = 3.
#       Así el código y la fórmula del paper se corresponden línea a línea,
#       sin argumentos algebraicos intermedios.
#
#   (b) LOG-LOSS. Se usa la entropía cruzada CATEGÓRICA: como H/D/A son
#       mutuamente excluyentes, solo importa la probabilidad asignada a la
#       categoría que efectivamente ocurrió. La entropía cruzada BINARIA
#       penalizaría también la probabilidad de las categorías que NO
#       ocurrieron y no corresponde a este problema.
#
#   (c) TEST DE SIGNIFICACIÓN DE LAS DIFERENCIAS (Sección 4).
#       Que un competidor tenga menor Brier Score que otro no dice si esa
#       ventaja es distinguible de cero, y esa es la pregunta directa de la
#       hipótesis de trabajo: ¿el mercado le gana significativamente al
#       modelo?
#
#       Los scores son promedios de 380 cantidades observadas partido a
#       partido, de modo que la diferencia entre dos competidores admite un
#       contraste pareado ordinario. Se implementa sobre las diferencias
#       individuales, con intervalo de confianza, contraste no paramétrico
#       de robustez y corrección por comparaciones múltiples.
#
#   (d) TRATAMIENTO DE FALTANTES EN LA SERIE DE CIERRE. Las métricas de la
#       serie de cierre se calculan sobre el subconjunto con cuota disponible
#       (datos_c / datos_largo_c) en lugar de sobre el total. Con cobertura
#       completa el resultado es idéntico; con faltantes, calcular sobre el
#       total devolvería NA sin avisar.
#
#   (e) EL LOG-LOSS DESAGREGADO NO ES UNA REGLA DE PUNTUACIÓN PROPIA
#       (Secciones 6 y 7). log_loss() filtra a ocurrio == 1. Aplicada a la
#       métrica GLOBAL eso es correcto: condiciona sobre el resultado
#       observado de cada partido y es exactamente la entropía cruzada
#       categórica. Pero aplicada DENTRO de un decil o de una categoría pasa
#       a puntuar únicamente las filas donde el evento ocurrió, y entonces la
#       esperanza se minimiza reportando f = 1: el predictor que INFLA
#       probabilidades gana, esté calibrado o no.
#
#       La demostración está en los propios datos. En el decil 1 la frecuencia
#       observada es 0,079; el mercado asigna 0,113 y el modelo 0,157. El
#       modelo está claramente PEOR calibrado ahí, el Brier Score —que sí es
#       propio— lo dice, y el log-loss condicional dice lo contrario
#       precisamente porque el modelo infla.
#
#       Por eso no se declara un ganador por log-loss desagregado. Las tablas
#       de log-loss por decil y por categoría se conservan como descriptivas
#       de la magnitud del error condicional, explícitamente marcadas como NO
#       aptas para declarar un ganador.
#
#       Regla general: condicionar sobre el PRONÓSTICO (los deciles) es
#       legítimo; condicionar sobre el RESULTADO (ocurrio == 1) rompe la
#       propiedad de la regla.
#
#   (f) REL DE MURPHY: PISO DE RUIDO (Sección 5). REL es un estimador sesgado
#       hacia arriba: con 10 bins de ~38 observaciones recoge el error de
#       muestreo de la frecuencia de cada bin como si fuera desajuste de
#       calibración. El valor esperado de REL para un predictor PERFECTAMENTE
#       calibrado es del orden de 0,005 por categoría y ~0,015 sumando las
#       tres, es decir del mismo orden que las REL observadas. Comparar las
#       REL de los tres competidores sería comparar ruido.
#
#       Por eso se estima el piso de ruido por simulación y REL no participa
#       del cuadro de ganadores. La conclusión de la Sección 5.1 —RES está muy
#       por encima del piso— queda sostenida cuantitativamente.
#
#   (g) EXCLUSIÓN DE LOS PARTIDOS DE LOS ASCENDIDOS (Sección 4.6).
#       Los equipos ascendidos tienen parámetros imputados vía
#       peso_descendidos (Script 04, Sección 3.1), un grado de libertad del
#       diseño sin justificación externa, y participan en ~10% de la
#       validación. Se recalcula todo sin ellos.
#
#   (h) ADVERTENCIAS HEREDADAS DECLARADAS AL PRINCIPIO (Sección 1.1).
#       Del Script 04: el cambio de régimen de la ventaja local debe
#       declararse ANTES de interpretar la calibración en la categoría H.
#       Del Script 07: la temporada de validación no reproduce el patrón de
#       calibración por deciles del período completo. Las dos se imprimen
#       antes de cualquier métrica, porque las dos cambian cómo hay que leer
#       todo lo que sigue.
#
#   (i) LECTURAS AUTOMÁTICAS. Los veredictos que imprime el script se redactan
#       para no afirmar más de lo que los datos sostienen. Dos casos:
#
#       (i1) Sección 1.1. Modelo y mercado sobreestiman P(H) casi en la misma
#       magnitud. El cambio de régimen de la ventaja local explica el NIVEL
#       del sesgo de ambos predictores y NO la BRECHA entre ellos: el mercado
#       tampoco anticipó la caída de 24/25. Atribuir el peor desempeño del
#       modelo en la categoría H al cambio de régimen sería una afirmación que
#       los datos no sostienen; esa diferencia tiene que venir de
#       discriminación entre partidos.
#
#       (i2) Sección 4.6. La desventaja del modelo no es "general": la brecha
#       está ausente en los partidos de ascendidos y concentrada en el resto.
#       Se reporta el error estándar de la submuestra —imprescindible, porque
#       son pocos partidos y no permiten afirmar una inversión— y la lectura
#       sustantiva: el modelo compite de igual a igual justamente donde el
#       mercado tiene menos información acumulada.
# -----------------------------------------------------------------------------

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. CARGAR DATOS
# -----------------------------------------------------------------------------
# Del Script 07 obtenemos, sobre la validación 2024/25:
#   - datos:       380 filas (formato ancho) con pi_* (pre-cierre),
#                  piC_* (cierre) y prob_* (modelo)
#   - datos_largo: 1.140 filas (380 partidos × 3 resultados) con:
#       pi_mercado         -> probabilidad implícita normalizada, pre-cierre
#       pi_mercado_cierre  -> probabilidad implícita normalizada, cierre
#       prob_modelo        -> probabilidad estimada por el modelo
#       ocurrio            -> 1 si ese resultado ocurrió, 0 si no
#
# La comparación queda entre TRES competidores:
#   1. Modelo Maher
#   2. Mercado — cuotas pre-cierre  (serie primaria)
#   3. Mercado — cuotas de cierre   (robustez)

mercado <- readRDS("resultados_mercado.rds")

datos_largo <- mercado$datos_largo
datos       <- mercado$datos

cat("Partidos:", nrow(datos), "| Filas en datos_largo:", nrow(datos_largo), "\n")

# --- Insumos heredados de los Scripts 01 y 04 --------------------------------
# Se cargan acá porque dos advertencias anotadas en scripts anteriores tienen
# que quedar declaradas ANTES de leer cualquier métrica de calibración, y
# porque el control de sensibilidad de la Sección 4.6 necesita saber qué
# partidos involucran a un equipo ascendido.
#
# Si alguno de los archivos no está, el script sigue corriendo y lo avisa: las
# secciones que dependen de ellos se saltean en lugar de romper la corrida.

equipos_nuevos <- if (file.exists("equipos_nuevos.rds")) {
  readRDS("equipos_nuevos.rds")
} else {
  cat("AVISO: no se encontró 'equipos_nuevos.rds' (Script 01).\n")
  character(0)
}

maher <- if (file.exists("resultados_maher.rds")) {
  readRDS("resultados_maher.rds")
} else {
  cat("AVISO: no se encontró 'resultados_maher.rds' (Script 04).\n")
  NULL
}

# Marca de los partidos que involucran a un equipo ascendido: son aquellos en
# los que el modelo opera con parámetros IMPUTADOS y no estimados.
datos <- datos |>
  mutate(ascendido = home_team %in% equipos_nuevos |
                     away_team %in% equipos_nuevos)

datos_largo <- datos_largo |>
  mutate(ascendido = home_team %in% equipos_nuevos |
                     away_team %in% equipos_nuevos)

if (length(equipos_nuevos) > 0) {
  cat("Equipos ascendidos (parámetros imputados):",
      paste(equipos_nuevos, collapse = ", "), "\n")
  cat(sprintf("Partidos de validación que los involucran: %d de %d (%.1f%%)\n",
              sum(datos$ascendido), nrow(datos),
              mean(datos$ascendido) * 100))
}

# --- Control de faltantes en la serie de cierre ---
# Si algún partido de 24/25 no tiene cuota de cierre, las métricas de esa
# serie se calculan sobre los partidos completos. Hay que declararlo, porque
# entonces las tres columnas de la tabla final no están sobre el mismo N.

n_sin_cierre <- sum(is.na(datos$piC_H) | is.na(datos$piC_D) | is.na(datos$piC_A))

if (n_sin_cierre > 0) {
  cat(sprintf(
    "\nATENCIÓN: %d partidos de 24/25 sin cuota de cierre (%.2f%%).\n",
    n_sin_cierre, n_sin_cierre / nrow(datos) * 100))
  cat("Las métricas de la serie de cierre se calculan sobre",
      nrow(datos) - n_sin_cierre, "partidos.\n")
  cat("DECLARAR ESTA DIFERENCIA DE N AL REPORTAR LA TABLA COMPARATIVA.\n")
} else {
  cat("\nTodos los partidos de 24/25 tienen cuota de cierre: las tres series\n")
  cat("se comparan sobre el mismo N =", nrow(datos), "partidos.\n")
}

# Subconjuntos completos para la serie de cierre
datos_c <- datos |>
  filter(!is.na(piC_H), !is.na(piC_D), !is.na(piC_A))

datos_largo_c <- datos_largo |>
  filter(!is.na(pi_mercado_cierre))

# Etiquetas usadas en todos los gráficos y tablas
etiquetas_series <- c(
  "Mercado (pre-cierre)" = "#e74c3c",
  "Mercado (cierre)"     = "#16a085",
  "Modelo Maher"         = "#2980b9"
)


# -----------------------------------------------------------------------------
# 1.1 ADVERTENCIAS HEREDADAS QUE CONDICIONAN LA LECTURA DE TODO EL SCRIPT
# -----------------------------------------------------------------------------
# Los Scripts 04 y 07 dejaron dos advertencias dirigidas explícitamente a este
# script. Se imprimen ACÁ, antes de cualquier métrica, porque las dos cambian
# cómo hay que leer los resultados que siguen. Dejarlas para el final sería
# invitar a que el lector saque la conclusión equivocada por el camino.

cat("\n", strrep("=", 74), "\n", sep = "")
cat("ADVERTENCIAS HEREDADAS DE LOS SCRIPTS 04 Y 07\n")
cat(strrep("=", 74), "\n", sep = "")

# --- (1) Cambio de régimen de la ventaja local (Script 04, Sección 5.1) ------
cat("\n--- (1) VENTAJA LOCAL: CAMBIO DE RÉGIMEN ENTRE MUESTRAS ---\n")

if (!is.null(maher) && !is.null(maher$ventaja_local_regimen)) {

  print(maher$ventaja_local_regimen |>
          mutate(across(where(is.double), \(x) round(x, 4))))

  cat(sprintf("\nBrecha del modelo respecto de la validación: %+.1f%%\n",
              maher$brecha_ventaja_pct))
  cat("\nEl modelo estima la ventaja local como una CONSTANTE sobre 21/22-23/24.\n")
  cat("En 24/25 el cociente observado fue menor, de modo que el modelo\n")
  cat("SOBREESTIMA P(H) de forma sistemática en los 380 partidos de validación.\n")

} else {
  cat("No se encontró 'resultados_maher.rds': no se puede reproducir la tabla\n")
  cat("de régimen. Correr el Script 04 antes que este.\n")
}

# Verificación directa sobre los datos de validación. No depende del Script 04:
# compara la probabilidad media que cada competidor asignó a cada categoría
# contra la frecuencia efectivamente observada.
sesgo_medio <- tibble(
  categoria      = c("H", "D", "A"),
  frecuencia_obs = c(mean(datos$ftr == "H"),
                     mean(datos$ftr == "D"),
                     mean(datos$ftr == "A")),
  media_mercado  = c(mean(datos$pi_H),   mean(datos$pi_D),   mean(datos$pi_A)),
  media_modelo   = c(mean(datos$prob_H), mean(datos$prob_D), mean(datos$prob_A))
) |>
  mutate(sesgo_mercado = media_mercado - frecuencia_obs,
         sesgo_modelo  = media_modelo  - frecuencia_obs)

cat("\n--- Probabilidad media asignada vs. frecuencia observada (24/25) ---\n")
print(sesgo_medio |> mutate(across(where(is.numeric), \(x) round(x, 4))))
cat("(sesgo positivo = el predictor asignó en promedio MÁS probabilidad a esa\n")
cat(" categoría de la que efectivamente se observó)\n")

# La pregunta no es solo si el modelo sobreestima P(H) —el Script 04 lo
# anticipaba— sino si lo hace MÁS que el mercado. Si los dos se equivocan en
# el nivel por la misma magnitud, el cambio de régimen explica el error de
# ambos y no la diferencia de desempeño entre ellos, que es lo que este
# script tiene que explicar.
sesgo_h_modelo  <- sesgo_medio$sesgo_modelo[1]
sesgo_h_mercado <- sesgo_medio$sesgo_mercado[1]
dif_sesgo_h     <- sesgo_h_modelo - sesgo_h_mercado

if (sesgo_h_modelo > 0 && sesgo_h_mercado > 0) {

  cat(sprintf(
    "\nLOS DOS sobreestiman P(H): el modelo en %+.4f y el mercado en %+.4f.\n",
    sesgo_h_modelo, sesgo_h_mercado))
  cat(sprintf("La diferencia entre ambos es de apenas %+.4f.\n", dif_sesgo_h))

  cat("\nEL CAMBIO DE RÉGIMEN LE PEGÓ A LOS DOS POR IGUAL. El mercado tampoco\n")
  cat("anticipó la caída de la ventaja local en 24/25. Es un dato con contenido\n")
  cat("propio sobre esa temporada y conviene reportarlo: el quiebre no era\n")
  cat("previsible con la información disponible, ni siquiera para quien fija\n")
  cat("los precios.\n")

  cat("\nCONSECUENCIA PARA LA REDACCIÓN — leer con atención:\n")
  cat("El cambio de régimen explica el NIVEL del sesgo de los dos predictores,\n")
  cat("pero NO explica la BRECHA entre ellos. NO escribir 'el peor desempeño\n")
  cat("del modelo en la categoría H se debe al cambio de régimen': si el\n")
  cat("mercado erró el nivel casi en la misma magnitud y aun así predice mejor\n")
  cat("esa categoría, la diferencia viene de DISCRIMINACIÓN entre partidos y no\n")
  cat("del nivel medio. Es exactamente lo que dice la Sección 5: la brecha es\n")
  cat("de resolución, no de calibración.\n")

} else if (sesgo_h_modelo > sesgo_h_mercado) {

  cat(sprintf(
    "\nSolo el modelo sobreestima P(H): %+.4f contra %+.4f del mercado.\n",
    sesgo_h_modelo, sesgo_h_mercado))
  cat("Acá el cambio de régimen SÍ diferencia a los dos predictores y puede\n")
  cat("usarse para explicar la brecha en la categoría H, no solo el nivel.\n")

} else {

  cat("\nLa sobreestimación de P(H) anticipada por el Script 04 NO aparece en el\n")
  cat("promedio de la validación. Revisar antes de usar ese argumento.\n")
}

cat("\nCONSECUENCIA PARA LA LECTURA DE ESTE SCRIPT:\n")
cat("  - El modelo tiene la ventaja local congelada en el promedio de\n")
cat("    21/22-23/24 y la validación operó por debajo de ese nivel. Es una\n")
cat("    limitación del DISEÑO TEMPORAL (parámetros estáticos), no un defecto\n")
cat("    de la especificación Poisson. Escribirlo así.\n")
cat("  - ALCANCE DEL ARGUMENTO: el cambio de régimen explica que el modelo\n")
cat("    haya errado el nivel de P(H). Si el mercado erró en magnitud similar\n")
cat("    (ver la comparación de arriba), NO explica además por qué el modelo\n")
cat("    predice peor esa categoría. Usar el argumento para el nivel, no para\n")
cat("    la brecha.\n")
cat("  - NO re-estimar la ventaja local con datos de 24/25: sería look-ahead\n")
cat("    bias y rompería el diseño fuera de muestra.\n")

# --- (2) El mercado de 24/25 no reproduce el patrón del período -------------
cat("\n--- (2) LA TEMPORADA DE VALIDACIÓN NO REPRODUCE EL PATRÓN DEL PERÍODO ---\n")

if (!is.null(mercado$correlaciones) && !is.null(mercado$datos_bins)) {

  cor_4temp <- mercado$correlaciones$correlacion[
    mercado$correlaciones$serie == "Pre-cierre"]
  cor_2425  <- cor(mercado$datos_bins$pi_mercado_med,
                   mercado$datos_bins$sesgo_mercado)

  cat(sprintf("Correlación pi-sesgo por decil, 4 temporadas (Script 07): %+.4f\n",
              cor_4temp))
  cat(sprintf("Correlación pi-sesgo por decil, solo 2024/25:             %+.4f\n",
              cor_2425))

  if (sign(cor_4temp) != sign(cor_2425)) {
    cat("\nEL SIGNO SE INVIERTE. Sobre las cuatro temporadas el patrón de\n")
    cat("deciles es compatible con el favourite-longshot bias clásico; sobre\n")
    cat("los 380 partidos de validación, no lo es.\n")
    cat("\nCONSECUENCIA: la comparación modelo vs. mercado de este script se\n")
    cat("hace sobre una temporada cuyo patrón de calibración NO es el del\n")
    cat("período. La ventaja del mercado que se mide acá describe 2024/25 y no\n")
    cat("se puede extrapolar hacia atrás ni presentarse como propiedad general\n")
    cat("del mercado. Declararlo junto con la tabla de resultados.\n")
  } else {
    cat("\nEl signo coincide con el del período completo: la validación\n")
    cat("reproduce cualitativamente el patrón general.\n")
  }

} else {
  cat("No están 'correlaciones' ni 'datos_bins' en resultados_mercado.rds:\n")
  cat("correr la versión actual del Script 07.\n")
}

# -----------------------------------------------------------------------------
# 2. FUNCIONES AUXILIARES
# -----------------------------------------------------------------------------
# Definimos las métricas como funciones para poder reutilizarlas
# en los distintos niveles de análisis (global, por resultado, por bin).

# r: cantidad de categorías de resultado (H, D, A).
r <- 3   # se conserva como constante documental; no interviene en el cálculo

# brier_global(): Brier Score clásico (Brier, 1950) — cálculo DIRECTO.
# Opera sobre datos (formato ancho, 380 filas), implementando literalmente
# la ecuación (2) de Brier (1950):
#
#   BS = (1/n) Σᵢ [ (f_H,i − o_H,i)² + (f_D,i − o_D,i)² + (f_A,i − o_A,i)² ]
#
# donde f_j,i es la probabilidad predicha para el resultado j en el partido i,
# y o_j,i = 1 si j ocurrió, 0 si no.
#
# Argumentos:
#   datos      → data frame con una fila por partido
#   ph, pd, pa → nombres de las columnas de probabilidad (H, D, A)
#
# Rango: [0, 2]. Menor es mejor.
brier_global <- function(datos, ph, pd, pa) {
  datos |>
    mutate(
      bs_i = (.data[[ph]] - as.integer(ftr == "H"))^2 +
             (.data[[pd]] - as.integer(ftr == "D"))^2 +
             (.data[[pa]] - as.integer(ftr == "A"))^2
    ) |>
    summarise(bs = mean(bs_i, na.rm = TRUE)) |>
    pull(bs)
}

# brier_score(): Brier Score binario — usado SOLO en secciones 6 y 7.
# Opera sobre vectores de datos_largo (formato largo).
# En esos contextos se aplica a un subconjunto de una sola categoría
# (p.ej. solo filas "H") o a filas agrupadas por decil — no a los 3
# resultados de un mismo partido. Es un Brier Score binario (rango [0,1])
# que sirve como diagnóstico de en qué resultado o rango se concentra el error.
# NO se usa para el cálculo global ni para la descomposición de Murphy.
brier_score <- function(prob, ocurrio) {
  mean((prob - ocurrio)^2, na.rm = TRUE)
}

# Log-loss / entropía cruzada categórica (pérdida informacional).
# Referencia: Good (1952) como fuente fundacional; Hvattum & Arntzen (2010)
# para su uso en predicción de fútbol.
#
# Fórmula: LL = -(1/n) * sum_i log(f_i, categoría ocurrida)
#
# Como H, D y A son mutuamente excluyentes, sólo hace falta evaluar la
# probabilidad que el predictor asignó a la categoría que realmente
# ocurrió en cada partido — a diferencia del Brier Score, acá NO se suma
# sobre las r categorías. Por eso la función filtra internamente a las
# filas donde ocurrio == 1 antes de promediar.
#
# Rango: [0, +∞). Menor es mejor.
# El clamp evita log(0) = -Inf cuando prob es exactamente 0.
log_loss <- function(prob, ocurrio) {
  prob_ocurrido <- prob[ocurrio == 1]
  prob_ocurrido <- prob_ocurrido[!is.na(prob_ocurrido)]
  prob_clamp    <- pmax(pmin(prob_ocurrido, 1 - 1e-10), 1e-10)
  -mean(log(prob_clamp))
}

# descomponer_brier(): descomposición de Murphy (1973), BS = REL − RES + UNC.
#
# La descomposición se aplica POR CATEGORÍA (H, D, A) en su versión binaria
# y luego se SUMA sobre las tres. Esa suma devuelve directamente la escala
# del Brier Score clásico [0, 2], de modo que los componentes son
# comparables con el BS global calculado por brier_global().
#
#   REL_j = (1/n_j) Σ_k n_jk (f_jk − o_jk)²      calibración   (menor mejor)
#   RES_j = (1/n_j) Σ_k n_jk (o_jk − ō_j)²       resolución    (mayor mejor)
#   UNC_j = ō_j (1 − ō_j)                        incertidumbre (fija)
#
# donde k indexa los bins dentro de la categoría j, f_jk es la probabilidad
# media predicha en el bin y o_jk la frecuencia observada en el bin.
#
# El binning se hace DENTRO de cada categoría: si se hiciera sobre el total,
# los bins mezclarían probabilidades de H, D y A, cuyas frecuencias base son
# distintas, y RES quedaría contaminada por esa diferencia.
#
# Los faltantes se descartan antes de binear: ntile() asigna NA a los
# valores faltantes y eso rompería el recuento por bin.
descomponer_brier <- function(datos_largo, col_prob, col_cat = "resultado",
                              col_obs = "ocurrio", n_bins = 10) {

  df <- datos_largo |>
    transmute(cat  = .data[[col_cat]],
              prob = .data[[col_prob]],
              obs  = .data[[col_obs]]) |>
    filter(!is.na(prob), !is.na(obs))

  # Frecuencia base y tamaño de CADA categoría (H, D, A) por separado
  bases <- df |>
    group_by(cat) |>
    summarise(n_cat = n(), base_cat = mean(obs), .groups = "drop")

  # Binning DENTRO de cada categoría y descomposición binaria por categoría
  comp <- df |>
    group_by(cat) |>
    mutate(bin = ntile(prob, n_bins)) |>
    group_by(cat, bin) |>
    summarise(n_bk = n(),
              f_bk = mean(prob),   # probabilidad media predicha en el bin
              o_bk = mean(obs),    # frecuencia real observada en el bin
              .groups = "drop") |>
    left_join(bases, by = "cat") |>
    group_by(cat) |>
    summarise(
      REL_j = sum(n_bk * (f_bk - o_bk)^2)     / first(n_cat),
      RES_j = sum(n_bk * (o_bk - base_cat)^2) / first(n_cat),
      UNC_j = first(base_cat) * (1 - first(base_cat)),
      .groups = "drop"
    )

  list(
    UNC = sum(comp$UNC_j),                       # = climatología
    REL = sum(comp$REL_j),
    RES = sum(comp$RES_j),
    BS_reconstruido = sum(comp$REL_j - comp$RES_j + comp$UNC_j),
    por_categoria = comp
  )
}

# scores_partido(): descompone las métricas globales en sus 380 sumandos.
#
# Tanto el Brier Score como el log-loss son PROMEDIOS de una cantidad
# definida partido a partido. Esta función devuelve esos sumandos, que son
# el insumo del test de significación de la Sección 4.
#
#   bs_i = Σ_j (f_ji − o_ji)²              contribución del partido i al BS
#   ll_i = −log(f_i, categoría ocurrida)   contribución del partido i al LL
#
# Por construcción: mean(bs_i) = BS global y mean(ll_i) = log-loss global.
# La Sección 4 verifica esa identidad antes de usar los vectores.
scores_partido <- function(d, ph, pd, pa) {

  o_H <- as.integer(d$ftr == "H")
  o_D <- as.integer(d$ftr == "D")
  o_A <- as.integer(d$ftr == "A")

  p_ocurrida <- case_when(
    d$ftr == "H" ~ d[[ph]],
    d$ftr == "D" ~ d[[pd]],
    d$ftr == "A" ~ d[[pa]]
  )

  tibble(
    bs_i = (d[[ph]] - o_H)^2 + (d[[pd]] - o_D)^2 + (d[[pa]] - o_A)^2,
    ll_i = -log(pmax(pmin(p_ocurrida, 1 - 1e-10), 1e-10))
  )
}


# -----------------------------------------------------------------------------
# 3. BRIER SCORE, BRIER SKILL SCORE Y LOG-LOSS — TRES COMPETIDORES
# -----------------------------------------------------------------------------
# Se comparan tres predictores sobre los mismos 380 partidos de validación:
#
#   1. Mercado (pre-cierre) — cuotas B365H/D/A, serie primaria
#   2. Mercado (cierre)     — cuotas B365CH/CD/CA, robustez
#   3. Modelo               — Maher puro (Poisson independiente)
#
# Las tres métricas se calculan igual para los tres:
#   - Brier Score: desde formato ancho con brier_global(), implementando la
#     ecuación (2) de Brier (1950). Rango [0, 2].
#   - Log-loss: entropía cruzada categórica sobre formato largo (la función
#     filtra internamente a ocurrio == 1).
#   - Brier Skill Score: mejora relativa sobre el predictor climatológico.

# --- Referencia climatológica del BSS ---
# Asigna a cada partido las frecuencias base observadas de 2024/25.
# Su Brier Score coincide con la incertidumbre (UNC) de la descomposición
# de Murphy: BS_ref = Σ_j base_j (1 − base_j).
#
# Es IDÉNTICA para los tres competidores. Si se recalculara por serie de
# cuotas, los BSS dejarían de ser comparables entre sí, que es justamente
# lo único para lo que sirve un skill score.
base_res <- prop.table(table(datos$ftr))
bs_ref   <- sum(base_res * (1 - base_res))

# --- Cálculo de las tres métricas para un competidor ---
# d       -> data frame en formato ancho sobre el que se calcula el BS
# dl      -> data frame en formato largo sobre el que se calcula el log-loss
# cH,cD,cA-> columnas de probabilidad en formato ancho
# c_larga -> columna de probabilidad en formato largo
#
# Cada competidor pasa su propio par (d, dl): así la serie de cierre se
# calcula sobre los partidos con cuota disponible y no devuelve NA.
metricas_de <- function(etiqueta, d, dl, cH, cD, cA, c_larga) {
  tibble(
    fuente      = etiqueta,
    n           = nrow(d),
    brier_score = brier_global(d, cH, cD, cA),
    log_loss    = log_loss(dl[[c_larga]], dl$ocurrio)
  )
}

metricas_globales <- bind_rows(
  metricas_de("Mercado (pre-cierre)", datos,   datos_largo,
              "pi_H",   "pi_D",   "pi_A",   "pi_mercado"),
  metricas_de("Mercado (cierre)",     datos_c, datos_largo_c,
              "piC_H",  "piC_D",  "piC_A",  "pi_mercado_cierre"),
  metricas_de("Modelo",               datos,   datos_largo,
              "prob_H", "prob_D", "prob_A", "prob_modelo")
) |>
  mutate(brier_skill = 1 - brier_score / bs_ref) |>
  select(fuente, n, brier_score, brier_skill, log_loss)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("MÉTRICAS GLOBALES — TRES COMPETIDORES (validación 2024/25)\n")
cat(strrep("=", 74), "\n", sep = "")
metricas_globales |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)
cat(sprintf("\nReferencia climatológica del BSS (idéntica para los tres): %.6f\n", bs_ref))

# --- Quién gana cada métrica ---
ganador <- function(valores, fuentes, direccion) {
  fuentes[if (direccion == "menor") which.min(valores) else which.max(valores)]
}

brecha <- function(valores, direccion) {
  ordenados <- sort(valores, decreasing = (direccion == "mayor"))
  abs(ordenados[1] - ordenados[2])   # distancia entre el 1° y el 2°
}

resumen_ganadores <- tibble(
  metrica   = c("Brier Score", "Brier Skill Score", "Log-loss"),
  direccion = c("menor", "mayor", "menor"),
  mejor     = c(
    ganador(metricas_globales$brier_score, metricas_globales$fuente, "menor"),
    ganador(metricas_globales$brier_skill, metricas_globales$fuente, "mayor"),
    ganador(metricas_globales$log_loss,    metricas_globales$fuente, "menor")
  ),
  ventaja_sobre_2do = c(
    brecha(metricas_globales$brier_score, "menor"),
    brecha(metricas_globales$brier_skill, "mayor"),
    brecha(metricas_globales$log_loss,    "menor")
  )
)

cat("\n--- Mejor competidor en cada métrica ---\n")
resumen_ganadores |>
  mutate(ventaja_sobre_2do = round(ventaja_sobre_2do, 6)) |>
  print(width = Inf)

cat("\nOJO: esta tabla dice QUIÉN tiene el mejor número, no si la ventaja es\n")
cat("distinguible de cero. Ese contraste es el de la Sección 4, y es el que\n")
cat("debe sostener cualquier afirmación del capítulo de resultados.\n")

# El Brier Skill Score es una transformación monótona decreciente del Brier
# Score (BSS = 1 − BS/BS_ref, con BS_ref constante), de modo que ambos DEBEN
# señalar al mismo ganador. Que coincidan no es un resultado: es un control
# de que la implementación está bien.
if (resumen_ganadores$mejor[1] != resumen_ganadores$mejor[2]) {
  warning("BS y BSS señalan ganadores distintos: imposible si BS_ref es único. Revisar.")
}

# El log-loss SÍ puede señalar otro ganador. Penaliza mucho más las
# probabilidades bajas asignadas al resultado que efectivamente ocurrió
# (crece sin cota cuando esa probabilidad tiende a cero), mientras que el
# Brier Score está acotado. Un predictor con pocos errores muy graves puede
# ganar en Brier y perder en log-loss.
if (resumen_ganadores$mejor[3] != resumen_ganadores$mejor[1]) {
  cat("\nATENCIÓN: el log-loss señala un ganador distinto al del Brier Score.\n")
  cat("No es un error. Indica que los competidores difieren en la GRAVEDAD de\n")
  cat("sus peores errores, no solo en el error promedio. Vale comentarlo.\n")
}

# --- Comparaciones de a pares ---
# Incluye explícitamente pre-cierre contra cierre, que es la pregunta de
# robustez: ¿cambia algo la conclusión al usar la otra serie de cuotas?
valor_de <- function(columna, fuente) {
  metricas_globales[[columna]][metricas_globales$fuente == fuente]
}

comparaciones_globales <- tibble(
  comparacion = c("Cierre − Pre-cierre",
                  "Modelo − Pre-cierre",
                  "Modelo − Cierre"),
  dif_brier = c(
    valor_de("brier_score", "Mercado (cierre)") - valor_de("brier_score", "Mercado (pre-cierre)"),
    valor_de("brier_score", "Modelo")           - valor_de("brier_score", "Mercado (pre-cierre)"),
    valor_de("brier_score", "Modelo")           - valor_de("brier_score", "Mercado (cierre)")
  ),
  dif_logloss = c(
    valor_de("log_loss", "Mercado (cierre)") - valor_de("log_loss", "Mercado (pre-cierre)"),
    valor_de("log_loss", "Modelo")           - valor_de("log_loss", "Mercado (pre-cierre)"),
    valor_de("log_loss", "Modelo")           - valor_de("log_loss", "Mercado (cierre)")
  )
) |>
  mutate(
    lectura_brier = if_else(dif_brier < 0,
                            "gana el primero", "gana el segundo")
  )

cat("\n--- Comparaciones de a pares (negativo = gana el primer término) ---\n")
comparaciones_globales |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

# --- Escalares para las Secciones 12 y 13 ---
bs_mercado  <- valor_de("brier_score", "Mercado (pre-cierre)")
bs_cierre   <- valor_de("brier_score", "Mercado (cierre)")
bs_modelo   <- valor_de("brier_score", "Modelo")
bss_mercado <- valor_de("brier_skill", "Mercado (pre-cierre)")
bss_cierre  <- valor_de("brier_skill", "Mercado (cierre)")
bss_modelo  <- valor_de("brier_skill", "Modelo")
ll_mercado  <- valor_de("log_loss", "Mercado (pre-cierre)")
ll_cierre   <- valor_de("log_loss", "Mercado (cierre)")
ll_modelo   <- valor_de("log_loss", "Modelo")


# =============================================================================
# 4. TEST DE SIGNIFICACIÓN DE LAS DIFERENCIAS ENTRE COMPETIDORES
# =============================================================================
# POR QUÉ HACE FALTA ESTA SECCIÓN.
#
# La Sección 3 informa que un competidor tiene menor Brier Score que otro,
# pero no si esa ventaja podría deberse al azar de una muestra de 380
# partidos. Sin este contraste, la afirmación central del trabajo —que el
# mercado supera al modelo, o que no lo supera— queda sin respaldo
# inferencial, y es lo primero que se pregunta al leer la tabla.
#
# EL ARGUMENTO.
#
# Tanto el Brier Score como el log-loss son PROMEDIOS MUESTRALES de una
# cantidad definida partido a partido:
#
#     BS = (1/n) Σᵢ bs_i        con bs_i = Σ_j (f_ji − o_ji)²
#     LL = (1/n) Σᵢ ll_i        con ll_i = −log(f_i, categoría ocurrida)
#
# Dos competidores predicen sobre LOS MISMOS partidos, de modo que sus
# scores están pareados. Definiendo la diferencia de pérdida
#
#     d_i = score_i(A) − score_i(B)
#
# la diferencia entre los scores globales es exactamente la media de los
# d_i, y contrastar H₀: E[d_i] = 0 es un test de una media sobre 380
# observaciones pareadas. Un valor negativo de la media indica que A pierde
# menos, es decir que A es el mejor predictor.
#
# ESTE ES EL TEST DE DIEBOLD Y MARIANO (1995) aplicado a corte transversal.
# En su formulación original la varianza de la media de d_i se estima con
# corrección HAC porque las predicciones son de una serie de tiempo con
# horizonte h > 1 y los errores se solapan. Acá las predicciones son de
# partidos distintos, no hay solapamiento y no hay razón para esperar
# autocorrelación en el orden de las fechas, de modo que la corrección HAC
# se reduce al error estándar ordinario y el estadístico DM coincide con el
# t de Student pareado. Conviene decirlo así en el texto: no es "un t-test
# casero", es el DM en el caso h = 1 sin autocorrelación, que es el uso
# estándar para comparar pronósticos probabilísticos fuera de muestra.
#
# TRES SALVAGUARDAS.
#
#   1. La distribución de d_i es asimétrica y de colas gruesas (los partidos
#      con resultado sorpresivo dominan). Con n = 380 el teorema central del
#      límite respalda el t, pero se agrega el test de rangos con signo de
#      Wilcoxon como contraste no paramétrico. Si ambos coinciden, la
#      conclusión no depende del supuesto de normalidad.
#
#   2. Se hacen tres comparaciones por métrica sobre los mismos datos. Se
#      reporta el p-valor corregido por el método de Holm, que controla la
#      tasa de error por familia sin suponer independencia entre contrastes.
#      Se informan ambos, crudo y corregido: la comparación "Modelo vs.
#      pre-cierre" está pre-especificada como contraste primario y no
#      necesita corrección; las otras dos son exploratorias y sí.
#
#   3. Se reporta el intervalo de confianza del 95%, no solo el p-valor. Con
#      380 partidos, un resultado no significativo puede convivir con un
#      intervalo lo bastante ancho como para no descartar diferencias
#      relevantes: en ese caso lo que corresponde escribir es que la muestra
#      no permite distinguir, no que los competidores son equivalentes.
# -----------------------------------------------------------------------------

# --- 4.1 SCORES PARTIDO A PARTIDO ---------------------------------------------
# Se construyen sobre el objeto datos completo (380 filas) para que los tres
# competidores queden alineados fila a fila. La serie de cierre tendrá NA en
# los partidos sin cuota; el pareo los descarta más abajo.

sc_precierre <- scores_partido(datos, "pi_H",   "pi_D",   "pi_A")
sc_cierre    <- scores_partido(datos, "piC_H",  "piC_D",  "piC_A")
sc_modelo    <- scores_partido(datos, "prob_H", "prob_D", "prob_A")

scores_individuales <- tibble(
  partido       = seq_len(nrow(datos)),
  ftr           = datos$ftr,
  bs_precierre  = sc_precierre$bs_i,
  bs_cierre     = sc_cierre$bs_i,
  bs_modelo     = sc_modelo$bs_i,
  ll_precierre  = sc_precierre$ll_i,
  ll_cierre     = sc_cierre$ll_i,
  ll_modelo     = sc_modelo$ll_i
)

# --- Control de consistencia con la Sección 3 ---
# El promedio de los scores individuales debe reproducir exactamente las
# métricas globales. Si no lo hace, los vectores no son los sumandos de esas
# métricas y el test estaría contrastando otra cosa.
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECCIÓN 4 — SIGNIFICACIÓN DE LAS DIFERENCIAS\n")
cat(strrep("=", 74), "\n", sep = "")

control_consistencia <- tibble(
  cantidad = c("BS pre-cierre", "BS cierre", "BS modelo",
               "LL pre-cierre", "LL cierre", "LL modelo"),
  desde_seccion_3 = c(bs_mercado, bs_cierre, bs_modelo,
                      ll_mercado, ll_cierre, ll_modelo),
  desde_scores_individuales = c(
    mean(scores_individuales$bs_precierre, na.rm = TRUE),
    mean(scores_individuales$bs_cierre,    na.rm = TRUE),
    mean(scores_individuales$bs_modelo,    na.rm = TRUE),
    mean(scores_individuales$ll_precierre, na.rm = TRUE),
    mean(scores_individuales$ll_cierre,    na.rm = TRUE),
    mean(scores_individuales$ll_modelo,    na.rm = TRUE)
  )
) |>
  mutate(
    discrepancia = abs(desde_seccion_3 - desde_scores_individuales),
    estado       = if_else(discrepancia < 1e-10, "OK", "REVISAR")
  )

cat("\n--- Control: el promedio de los scores individuales reproduce la métrica global ---\n")
control_consistencia |>
  mutate(across(where(is.numeric), ~ signif(., 8))) |>
  print(width = Inf)

if (any(control_consistencia$estado == "REVISAR")) {
  warning("Los scores individuales no promedian a las métricas globales. ",
          "Revisar scores_partido() antes de interpretar los tests.")
}

# --- 4.2 FUNCIÓN DE CONTRASTE PAREADO -----------------------------------------
# H0: E[score(A) − score(B)] = 0
#
# Devuelve la diferencia media, su error estándar, el intervalo de confianza
# del 95%, el estadístico y el p-valor del t pareado, el p-valor de Wilcoxon
# y una medida de tamaño de efecto (d de Cohen para muestras pareadas).
#
# La proporción de partidos en que A pierde menos que B se reporta como
# lectura descriptiva: un competidor puede ganar en promedio y sin embargo
# ser peor en la mayoría de los partidos, si su ventaja se concentra en
# pocos casos de error grave.

#
# NOTA DE IMPLEMENTACIÓN. Los nombres de los competidores se pasan como
# argumentos separados (nombre_A, nombre_B) en lugar de extraerse de la
# etiqueta con expresiones regulares. Parsear la etiqueta obligaría a buscar
# el signo menos, que acá es el carácter tipográfico U+2212 y no el guion
# ASCII: cualquier apertura del archivo con un encoding distinto rompería el
# filtro en silencio y devolvería cero filas. Con argumentos explícitos el
# código no depende del encoding.

comparar_scores <- function(s_A, s_B, nombre_A, nombre_B, metrica) {

  dif  <- s_A - s_B
  usar <- !is.na(dif)
  dif  <- dif[usar]

  tt <- t.test(dif)
  wt <- suppressWarnings(wilcox.test(dif, exact = FALSE))

  tibble(
    metrica      = metrica,
    competidor_A = nombre_A,
    competidor_B = nombre_B,
    comparacion  = paste(nombre_A, "-", nombre_B),
    n_pares      = length(dif),
    dif_media    = mean(dif),
    ee           = sd(dif) / sqrt(length(dif)),
    ic95_inf     = unname(tt$conf.int[1]),
    ic95_sup     = unname(tt$conf.int[2]),
    t            = unname(tt$statistic),
    gl           = unname(tt$parameter),
    p_t          = tt$p.value,
    p_wilcoxon   = wt$p.value,
    d_cohen      = mean(dif) / sd(dif),
    prop_gana_A  = mean(dif < 0)
  )
}

# --- 4.3 LOS SEIS CONTRASTES --------------------------------------------------
# Tres comparaciones × dos métricas. El orden de los términos es siempre
# "primero − segundo", de modo que una diferencia NEGATIVA indica que el
# PRIMER competidor pierde menos, es decir que es el mejor.

tests_brier <- bind_rows(
  comparar_scores(scores_individuales$bs_modelo,
                  scores_individuales$bs_precierre,
                  "Modelo", "Pre-cierre", "Brier Score"),
  comparar_scores(scores_individuales$bs_modelo,
                  scores_individuales$bs_cierre,
                  "Modelo", "Cierre",     "Brier Score"),
  comparar_scores(scores_individuales$bs_cierre,
                  scores_individuales$bs_precierre,
                  "Cierre", "Pre-cierre", "Brier Score")
)

tests_logloss <- bind_rows(
  comparar_scores(scores_individuales$ll_modelo,
                  scores_individuales$ll_precierre,
                  "Modelo", "Pre-cierre", "Log-loss"),
  comparar_scores(scores_individuales$ll_modelo,
                  scores_individuales$ll_cierre,
                  "Modelo", "Cierre",     "Log-loss"),
  comparar_scores(scores_individuales$ll_cierre,
                  scores_individuales$ll_precierre,
                  "Cierre", "Pre-cierre", "Log-loss")
)

# Corrección de Holm DENTRO de cada métrica (familia = las 3 comparaciones)
tests_significacion <- bind_rows(tests_brier, tests_logloss) |>
  mutate(p_holm = p.adjust(p_t, method = "holm"), .by = metrica) |>
  mutate(
    significativo_05 = p_t < 0.05,
    mejor = case_when(
      p_t >= 0.05    ~ "indistinguibles",
      dif_media < 0  ~ competidor_A,
      TRUE           ~ competidor_B
    )
  ) |>
  select(metrica, competidor_A, competidor_B, comparacion,
         n_pares, dif_media, ee,
         ic95_inf, ic95_sup, t, gl, p_t, p_holm, p_wilcoxon,
         d_cohen, prop_gana_A, significativo_05, mejor)

cat("\n--- Contrastes pareados sobre los scores partido a partido ---\n")
cat("H0: la diferencia media de pérdida es cero\n")
cat("Diferencia NEGATIVA = el primer competidor pierde menos (es mejor)\n\n")

tests_significacion |>
  select(metrica, comparacion, n_pares, dif_media, ee, ic95_inf, ic95_sup,
         p_t, p_holm, p_wilcoxon) |>
  mutate(across(where(is.numeric), ~ signif(., 4))) |>
  print(n = Inf, width = Inf)

cat("\n--- Tamaño de efecto y lectura descriptiva ---\n")
tests_significacion |>
  select(metrica, comparacion, dif_media, d_cohen, prop_gana_A, mejor) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(n = Inf, width = Inf)

cat("\nprop_gana_A = proporción de partidos en que el PRIMER competidor tiene\n")
cat("menor pérdida. Si está cerca de 0,5 mientras la diferencia media es\n")
cat("significativa, la ventaja no es parejo partido a partido sino que se\n")
cat("concentra en un subconjunto de casos. Vale comentarlo en el texto.\n")

# --- 4.4 COHERENCIA ENTRE EL t Y EL WILCOXON ----------------------------------
# Se comparan las CONCLUSIONES, no los p-valores. Que uno dé 0,043 y el otro
# 0,061 no es una contradicción: es la misma evidencia con dos estimadores de
# precisión distintos. Solo importa si uno rechaza rotundamente y el otro no
# se acerca.

cat("\n--- Coherencia entre el contraste paramétrico y el no paramétrico ---\n")

coherencia <- tests_significacion |>
  transmute(
    metrica, comparacion,
    p_t        = round(p_t, 4),
    p_wilcoxon = round(p_wilcoxon, 4),
    acuerdo = case_when(
      (p_t < 0.05) == (p_wilcoxon < 0.05)          ~ "coinciden",
      abs(p_t - p_wilcoxon) < 0.05                 ~ "discrepan al margen",
      TRUE                                          ~ "DISCREPAN"
    )
  )

print(coherencia, n = Inf, width = Inf)

if (any(coherencia$acuerdo == "DISCREPAN")) {
  cat("\nATENCIÓN: hay al menos un contraste donde el t y el Wilcoxon difieren\n")
  cat("de forma sustantiva. Causa probable: asimetría fuerte de las\n")
  cat("diferencias, con unos pocos partidos de error extremo que mueven la\n")
  cat("media pero no la mediana. En ese caso el Wilcoxon es la referencia y\n")
  cat("hay que reportar ambos, explicando por qué difieren.\n")
} else {
  cat("\nLas dos vías coinciden: la conclusión no depende del supuesto de\n")
  cat("normalidad de las diferencias.\n")
}

# --- 4.5 CONCLUSIÓN DEL CONTRASTE PRIMARIO ------------------------------------
# El contraste PRE-ESPECIFICADO del trabajo es Modelo vs. Mercado pre-cierre
# en Brier Score. Los otros cinco son robustez (serie de cierre) y métrica
# alternativa (log-loss).

primario <- tests_significacion |>
  filter(metrica == "Brier Score",
         competidor_A == "Modelo", competidor_B == "Pre-cierre")

stopifnot("No se encontró el contraste primario" = nrow(primario) == 1)

cat("\n", strrep("-", 74), "\n", sep = "")
cat("CONTRASTE PRIMARIO: Modelo vs. Mercado (pre-cierre) en Brier Score\n")
cat(strrep("-", 74), "\n", sep = "")

cat(sprintf("\nDiferencia media (modelo − mercado): %+.6f\n", primario$dif_media))
cat(sprintf("IC 95%%: [%+.6f ; %+.6f]\n", primario$ic95_inf, primario$ic95_sup))
cat(sprintf("t(%d) = %.4f | p = %.4f | Wilcoxon p = %.4f | d = %.3f\n",
            primario$gl, primario$t, primario$p_t,
            primario$p_wilcoxon, primario$d_cohen))

if (primario$p_t < 0.05 && primario$dif_media > 0) {
  cat("\nEl mercado supera al modelo de forma estadísticamente significativa.\n")
  cat("REDACTAR: la ventaja del mercado no es atribuible al azar muestral.\n")
  cat("La Sección 5 dirá de dónde viene esa ventaja (calibración o resolución).\n")
} else if (primario$p_t < 0.05 && primario$dif_media < 0) {
  cat("\nEl MODELO supera al mercado de forma estadísticamente significativa.\n")
  cat("REDACTAR con cautela y verificar antes que no haya fuga de información:\n")
  cat("el modelo se estimó con 21/22–23/24 y se evalúa en 24/25, de modo que\n")
  cat("no debería tener acceso a información que el mercado no tuviera.\n")
} else {
  cat("\nLa diferencia NO alcanza significación al 5%.\n")
  cat("REDACTAR: la muestra de 380 partidos no permite distinguir el\n")
  cat("desempeño del modelo del desempeño del mercado en Brier Score.\n")
  cat("NO escribir que son equivalentes: mirar el ancho del intervalo.\n")
}

# --- Lectura por precisión del intervalo ---
# Mismo criterio que en la Sección 7.4 del Script 07: qué valores descarta el
# intervalo importa más que si el p-valor cruzó un umbral. Se expresa la
# amplitud como porcentaje del Brier Score del mercado, que es la escala en
# que el lector puede juzgar si una diferencia es relevante.
amplitud_rel <- (primario$ic95_sup - primario$ic95_inf) / bs_mercado * 100

cat(sprintf("\nAmplitud del IC 95%%: %.6f, equivalente al %.2f%% del Brier Score\n",
            primario$ic95_sup - primario$ic95_inf, amplitud_rel))
cat(sprintf("del mercado (%.6f).\n", bs_mercado))
cat("El intervalo descarta toda diferencia fuera de ese rango. Si es ancho,\n")
cat("la no significación indica falta de precisión, no equivalencia.\n")

# --- Robustez: ¿cambia la conclusión con la serie de cierre? ---
robustez_cierre <- tests_significacion |>
  filter(metrica == "Brier Score",
         competidor_A == "Modelo", competidor_B == "Cierre")

cat("\n--- Robustez: el mismo contraste con la serie de cierre ---\n")
cat(sprintf("Diferencia media: %+.6f | p = %.4f | n = %d pares\n",
            robustez_cierre$dif_media, robustez_cierre$p_t,
            robustez_cierre$n_pares))

if ((primario$p_t < 0.05) == (robustez_cierre$p_t < 0.05) &&
    sign(primario$dif_media) == sign(robustez_cierre$dif_media)) {
  cat("La conclusión NO depende de qué serie de cuotas se use: mismo signo y\n")
  cat("misma decisión en ambas. Reportarlo como control de robustez superado.\n")
} else {
  cat("ATENCIÓN: la conclusión CAMBIA según la serie de cuotas. Hay que\n")
  cat("reportarlo explícitamente y discutirlo: la serie de cierre incorpora\n")
  cat("información que la pre-cierre todavía no tenía, de modo que una mejora\n")
  cat("del mercado al cerrar es sustantiva, no un artefacto.\n")
}

# --- Concordancia con la métrica alternativa ---
primario_ll <- tests_significacion |>
  filter(metrica == "Log-loss",
         competidor_A == "Modelo", competidor_B == "Pre-cierre")

cat("\n--- Concordancia con el log-loss ---\n")
cat(sprintf("Brier:    diferencia %+.6f | p = %.4f\n",
            primario$dif_media, primario$p_t))
cat(sprintf("Log-loss: diferencia %+.6f | p = %.4f\n",
            primario_ll$dif_media, primario_ll$p_t))

if (sign(primario$dif_media) == sign(primario_ll$dif_media)) {
  cat("Las dos reglas de puntuación propias apuntan en la misma dirección.\n")
  cat("Es lo esperable y refuerza la conclusión.\n")
} else {
  cat("Las dos métricas apuntan en direcciones OPUESTAS. No es un error: el\n")
  cat("log-loss no está acotado y castiga con mucha más fuerza los errores\n")
  cat("graves. Indica que los competidores difieren en el perfil de sus\n")
  cat("peores fallos. Reportarlo y explicarlo.\n")
}


# -----------------------------------------------------------------------------
# 4.6 ROBUSTEZ: EXCLUSIÓN DE LOS PARTIDOS DE LOS EQUIPOS ASCENDIDOS
# -----------------------------------------------------------------------------
# Pendiente anotada en el Script 04, Sección 3.1.
#
# Los equipos ascendidos a 24/25 no tienen historial en la muestra de
# entrenamiento: sus parámetros de ataque y defensa no se estiman, se IMPUTAN
# mediante peso_descendidos, que es un grado de libertad del diseño sin
# justificación externa. Como intervienen en aproximadamente el 10% de los
# partidos de validación, el modelo compite en desventaja en esa fracción por
# una decisión arbitraria del diseño y no por su especificación. El mercado,
# en cambio, cotiza esos partidos con toda la información disponible.
#
# EL CONTROL. Recalcular las métricas y el contraste primario excluyendo esos
# partidos. Si la conclusión se mantiene, la elección del peso deja de ser
# relevante para los resultados y se declara resuelta. Si cambia, es una
# limitación que hay que escribir.
#
# ADVERTENCIA. Esto es un análisis de SENSIBILIDAD, no la especificación
# principal. Los partidos de los ascendidos son parte legítima de la temporada
# y excluirlos del resultado principal sería seleccionar la muestra en función
# de lo que conviene al modelo. El resultado que se reporta es el de la
# muestra completa; este es el control que lo acompaña.

cat("\n", strrep("=", 74), "\n", sep = "")
cat("4.6 — SENSIBILIDAD: EXCLUYENDO LOS PARTIDOS DE LOS ASCENDIDOS\n")
cat(strrep("=", 74), "\n", sep = "")

if (length(equipos_nuevos) == 0 || sum(datos$ascendido) == 0) {

  cat("\nNo hay partidos con equipos ascendidos: nada que excluir.\n")
  metricas_sin_asc   <- NULL
  primario_sin_asc   <- NULL
  comparacion_asc    <- NULL

} else {

  idx_sin_asc         <- !datos$ascendido
  datos_sin_asc       <- datos       |> filter(!ascendido)
  datos_largo_sin_asc <- datos_largo |> filter(!ascendido)

  cat(sprintf("\nPartidos excluidos: %d | partidos restantes: %d\n",
              sum(datos$ascendido), nrow(datos_sin_asc)))

  # --- Métricas globales sobre la submuestra ---
  metricas_sin_asc <- bind_rows(
    metricas_de("Mercado (pre-cierre)", datos_sin_asc, datos_largo_sin_asc,
                "pi_H", "pi_D", "pi_A", "pi_mercado"),
    metricas_de("Mercado (cierre)",
                datos_sin_asc       |> filter(!is.na(piC_H)),
                datos_largo_sin_asc |> filter(!is.na(pi_mercado_cierre)),
                "piC_H", "piC_D", "piC_A", "pi_mercado_cierre"),
    metricas_de("Modelo", datos_sin_asc, datos_largo_sin_asc,
                "prob_H", "prob_D", "prob_A", "prob_modelo")
  )

  comparacion_metricas <- metricas_globales |>
    select(fuente, bs_completo = brier_score, ll_completo = log_loss) |>
    left_join(
      metricas_sin_asc |> select(fuente, bs_sin_asc = brier_score,
                                 ll_sin_asc = log_loss),
      by = "fuente"
    ) |>
    mutate(dif_bs = bs_sin_asc - bs_completo,
           dif_ll = ll_sin_asc - ll_completo)

  cat("\n--- Métricas globales con y sin los partidos de ascendidos ---\n")
  comparacion_metricas |>
    mutate(across(where(is.numeric), ~ round(., 6))) |>
    print(width = Inf)

  # --- ¿Le va peor al modelo justamente en esos partidos? ---
  # Es la pregunta concreta que motivaba el control: si la brecha entre modelo
  # y mercado es mucho mayor en los partidos de ascendidos, la imputación está
  # contaminando el resultado global.
  comparacion_asc <- tibble(
    grupo      = c("Con ascendido", "Sin ascendido"),
    partidos   = c(sum(datos$ascendido), sum(!datos$ascendido)),
    bs_mercado = c(mean(scores_individuales$bs_precierre[ datos$ascendido]),
                   mean(scores_individuales$bs_precierre[!datos$ascendido])),
    bs_modelo  = c(mean(scores_individuales$bs_modelo[ datos$ascendido]),
                   mean(scores_individuales$bs_modelo[!datos$ascendido]))
  ) |>
    mutate(brecha_modelo_mercado = bs_modelo - bs_mercado)

  cat("\n--- Brecha modelo-mercado dentro y fuera de los partidos de ascendidos ---\n")
  comparacion_asc |>
    mutate(across(where(is.numeric), ~ round(., 6))) |>
    print(width = Inf)

  # --- Contraste primario sobre la submuestra ---
  primario_sin_asc <- comparar_scores(
    scores_individuales$bs_modelo[idx_sin_asc],
    scores_individuales$bs_precierre[idx_sin_asc],
    "Modelo", "Pre-cierre", "Brier Score")

  primario_ll_sin_asc <- comparar_scores(
    scores_individuales$ll_modelo[idx_sin_asc],
    scores_individuales$ll_precierre[idx_sin_asc],
    "Modelo", "Pre-cierre", "Log-loss")

  cat("\n--- Contraste primario: muestra completa vs. sin ascendidos ---\n")
  bind_rows(
    primario          |> mutate(muestra = "Completa"),
    primario_sin_asc  |> mutate(muestra = "Sin ascendidos"),
    primario_ll       |> mutate(muestra = "Completa"),
    primario_ll_sin_asc |> mutate(muestra = "Sin ascendidos")
  ) |>
    select(metrica, muestra, n_pares, dif_media, ic95_inf, ic95_sup, p_t) |>
    mutate(across(where(is.numeric), ~ signif(., 4))) |>
    print(n = Inf, width = Inf)

  # --- Veredicto ---
  cambia_decision <- (primario$p_t < 0.05) != (primario_sin_asc$p_t < 0.05)
  cambia_signo    <- sign(primario$dif_media) != sign(primario_sin_asc$dif_media)

  cat("\n")
  if (!cambia_decision && !cambia_signo) {
    cat("La conclusión NO cambia al excluir los partidos de los ascendidos:\n")
    cat("mismo signo y misma decisión al 5%. La elección de peso_descendidos\n")
    cat("del Script 04 no afecta el resultado y queda DECLARADA COMO RESUELTA:\n")
    cat("ya no es una limitación que haya que arrastrar al capítulo de\n")
    cat("conclusiones, es un grado de libertad cuyo efecto se midió y es nulo.\n")
    cat(sprintf("La ventaja del mercado pasa de %+.6f a %+.6f al excluirlos.\n",
                primario$dif_media, primario_sin_asc$dif_media))
  } else {
    cat("ATENCIÓN: la conclusión CAMBIA al excluir los partidos de los\n")
    cat("ascendidos. La imputación de parámetros vía peso_descendidos está\n")
    cat("influyendo sobre el resultado principal. Reportar las dos versiones y\n")
    cat("declararlo como limitación explícita del diseño: el modelo no puede\n")
    cat("evaluarse en igualdad de condiciones sobre equipos sin historial.\n")
  }

  if (!is.null(comparacion_asc)) {

    brecha_con <- comparacion_asc$brecha_modelo_mercado[1]
    brecha_sin <- comparacion_asc$brecha_modelo_mercado[2]

    # Error estándar de la brecha DENTRO de la submuestra de ascendidos. Sin
    # esto no se puede decir nada sobre ella: son pocos partidos y cualquier
    # lectura direccional tiene que declarar su propia imprecisión.
    dif_asc <- scores_individuales$bs_modelo[datos$ascendido] -
               scores_individuales$bs_precierre[datos$ascendido]
    ee_asc  <- sd(dif_asc) / sqrt(length(dif_asc))

    cat(sprintf("\nBrecha modelo-mercado: %+.6f en los %d partidos con ascendido\n",
                brecha_con, sum(datos$ascendido)))
    cat(sprintf("   (EE = %.6f -> %.2f errores estándar de cero)\n",
                ee_asc, abs(brecha_con) / ee_asc))
    cat(sprintf("                       %+.6f en los %d restantes.\n",
                brecha_sin, sum(!datos$ascendido)))

    if (brecha_sin > 0 && brecha_con > 2 * brecha_sin) {

      cat("\nLa desventaja del modelo se CONCENTRA en los partidos de equipos con\n")
      cat("parámetros imputados. Hay que decirlo y tratarlo como limitación del\n")
      cat("diseño: el modelo no se evalúa en igualdad de condiciones ahí.\n")

    } else if (brecha_sin > 0 && brecha_con < brecha_sin / 2) {

      cat("\nLa desventaja del modelo NO está repartida: está AUSENTE en esos\n")
      cat("partidos y concentrada en el resto. No corresponde escribir que es\n")
      cat("'general'.\n")

      if (abs(brecha_con) < 2 * ee_asc) {
        cat("\nCAUTELA: con esa cantidad de partidos, la brecha de la submuestra\n")
        cat("está a menos de dos errores estándar de cero. NO se puede afirmar\n")
        cat("que se invierta ni que el modelo sea mejor ahí. Lo que SÍ queda\n")
        cat("establecido —y es lo que importa para el pendiente del Script 04—\n")
        cat("es que la desventaja global no proviene de esos partidos.\n")
      }

      cat("\nLECTURA SUSTANTIVA (va al texto): el modelo compite de igual a igual\n")
      cat("justamente donde el mercado tiene menos información acumulada —un\n")
      cat("equipo recién ascendido, sin historial en la categoría— y pierde\n")
      cat("donde el mercado dispone de tres temporadas de antecedentes. Es\n")
      cat("consistente con que la ventaja del mercado sea de RESOLUCIÓN y de\n")
      cat("origen informacional, y no de calibración.\n")

    } else {

      cat("\nLa desventaja del modelo está repartida de forma pareja entre los dos\n")
      cat("grupos: es general y atribuible al modelo, no a la imputación.\n")
    }
  }
}


# -----------------------------------------------------------------------------
# 5. DESCOMPOSICIÓN DEL BRIER SCORE (Murphy, 1973)
# -----------------------------------------------------------------------------
#   BS = REL − RES + UNC
#
# INCERTIDUMBRE (UNC): dificultad intrínseca del problema. Depende solo de
#   las frecuencias observadas, no del predictor. Idéntica para los tres.
#
# CALIBRACIÓN (REL): ¿las probabilidades predichas coinciden con las
#   frecuencias reales dentro de cada bin? Mide sesgo. MENOR ES MEJOR
#   (un predictor perfectamente calibrado tiene REL = 0).
#
# RESOLUCIÓN (RES): ¿el predictor se aparta de la frecuencia base cuando
#   corresponde? Mide capacidad de discriminación. MAYOR ES MEJOR.
#
# La descomposición se calcula POR CATEGORÍA (H, D, A) y se suma sobre las
# tres, lo que devuelve directamente la escala del BS clásico [0, 2].
#
# La Sección 4 dijo SI hay diferencia; esta dice DE DÓNDE viene.

descomp_mercado <- descomponer_brier(datos_largo,   "pi_mercado")
descomp_cierre  <- descomponer_brier(datos_largo_c, "pi_mercado_cierre")
descomp_modelo  <- descomponer_brier(datos_largo,   "prob_modelo")

descomposicion <- tibble(
  fuente        = c("Mercado (pre-cierre)", "Mercado (cierre)", "Modelo"),
  calibracion   = c(descomp_mercado$REL, descomp_cierre$REL, descomp_modelo$REL),
  resolucion    = c(descomp_mercado$RES, descomp_cierre$RES, descomp_modelo$RES),
  incertidumbre = c(descomp_mercado$UNC, descomp_cierre$UNC, descomp_modelo$UNC),
  bs_reconstruido = c(descomp_mercado$BS_reconstruido,
                      descomp_cierre$BS_reconstruido,
                      descomp_modelo$BS_reconstruido),
  bs_calculado    = c(bs_mercado, bs_cierre, bs_modelo)
)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("DESCOMPOSICIÓN DEL BRIER SCORE (Murphy, 1973) — TRES COMPETIDORES\n")
cat(strrep("=", 74), "\n", sep = "")
descomposicion |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

# --- Piso de ruido de REL -----------------------------------------------------
# REL es un estimador SESGADO HACIA ARRIBA. Con 10 bins de ~38 observaciones,
# la frecuencia observada de cada bin tiene error de muestreo, y REL recoge ese
# error como si fuera desajuste de calibración. El valor esperado de REL para
# un predictor PERFECTAMENTE CALIBRADO no es cero: es aproximadamente
# (1/n) Σ_k p_k(1 − p_k), del orden de 0,005 por categoría a este tamaño de
# muestra y ~0,015 sumando las tres — es decir, del mismo orden que las REL
# efectivamente observadas.
#
# Consecuencia: comparar las REL de los tres competidores entre sí es comparar
# ruido de binning. Por eso REL sale del cuadro de "ganadores" y se reemplaza
# por un contraste contra su propio piso.
#
# LA SIMULACIÓN. Se generan resultados artificiales y ~ Bernoulli(f_i) usando
# las probabilidades REALES de cada competidor: por construcción ese predictor
# está perfectamente calibrado. Se recalcula REL con el mismo binning. La
# distribución resultante es el piso de ruido de este diseño.

rel_simulada <- function(dl, col_prob, n_sim = 200, semilla = 904158) {
  set.seed(semilla)
  base <- dl |>
    transmute(cat = resultado, prob = .data[[col_prob]]) |>
    filter(!is.na(prob))
  vapply(seq_len(n_sim), function(s) {
    sim <- base |> mutate(ocurrio_sim = rbinom(n(), 1, prob))
    descomponer_brier(sim, "prob", col_cat = "cat",
                      col_obs = "ocurrio_sim")$REL
  }, numeric(1))
}

cat("\n--- Piso de ruido de REL (simulación bajo calibración perfecta) ---\n")
cat("Puede tardar unos segundos...\n")

sim_rel_pre <- rel_simulada(datos_largo,   "pi_mercado")
sim_rel_cie <- rel_simulada(datos_largo_c, "pi_mercado_cierre")
sim_rel_mod <- rel_simulada(datos_largo,   "prob_modelo")

piso_rel <- tibble(
  fuente        = c("Mercado (pre-cierre)", "Mercado (cierre)", "Modelo"),
  REL_observada = c(descomp_mercado$REL, descomp_cierre$REL, descomp_modelo$REL),
  REL_esperada_si_calibrado = c(mean(sim_rel_pre), mean(sim_rel_cie),
                                mean(sim_rel_mod)),
  p95_ruido = c(unname(quantile(sim_rel_pre, 0.95)),
                unname(quantile(sim_rel_cie, 0.95)),
                unname(quantile(sim_rel_mod, 0.95)))
) |>
  mutate(excede_ruido = REL_observada > p95_ruido)

piso_rel |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

if (any(piso_rel$excede_ruido)) {
  cat("\nAl menos un competidor tiene REL por encima del percentil 95 del ruido:\n")
  cat("ahí sí hay desajuste de calibración detectable y se puede comentar.\n")
} else {
  cat("\nNINGÚN competidor supera el percentil 95 del piso de ruido. A este\n")
  cat("tamaño de muestra y con este binning, los tres son indistinguibles de\n")
  cat("un predictor perfectamente calibrado.\n")
  cat("\nESTE ES UN RESULTADO REPORTABLE, no una limitación: la diferencia de\n")
  cat("Brier Score entre modelo y mercado NO viene de calibración. Viene\n")
  cat("enteramente de resolución. Y refuerza la conclusión de la Sección 5.1\n")
  cat("en lugar de debilitarla.\n")
}

cat("\nPOR ESO REL NO ENTRA EN EL CUADRO DE GANADORES QUE SIGUE: ordenar a los\n")
cat("tres competidores por REL sería ordenar ruido de muestreo.\n")

# --- Quién gana el componente interpretable ---
# Solo RES. Sus diferencias (del orden de 0,02-0,04) están muy por encima del
# piso de ruido de la descomposición, de modo que el orden sí es informativo.
resumen_descomposicion <- tibble(
  componente = "Resolución (RES)",
  direccion  = "mayor",
  mejor      = ganador(descomposicion$resolucion, descomposicion$fuente, "mayor"),
  ventaja_sobre_2do = brecha(descomposicion$resolucion, "mayor")
)

cat("\n--- Mejor competidor en resolución ---\n")
resumen_descomposicion |>
  mutate(ventaja_sobre_2do = round(ventaja_sobre_2do, 6)) |>
  print(width = Inf)

cat("\nLECTURA. La descomposición dice DE DÓNDE viene la diferencia en el\n")
cat("Brier Score global. Que un predictor gane en REL y pierda en RES\n")
cat("significa que está mejor calibrado pero discrimina menos: sus\n")
cat("probabilidades son más honestas, pero se juegan menos lejos de la\n")
cat("frecuencia base. Es la distinción sustantiva del capítulo.\n")

cat("\nADVERTENCIA PARA LA REDACCIÓN. Los componentes REL y RES no llevan\n")
cat("error estándar en esta implementación: dependen del binning y su\n")
cat("distribución muestral no es la de una media simple. La Sección 4\n")
cat("contrasta el BS global, que sí es una media; la descomposición se\n")
cat("reporta como lectura DESCRIPTIVA del origen de la diferencia. No\n")
cat("escribir 'la ventaja en resolución es significativa' sin un contraste\n")
cat("propio (requeriría bootstrap sobre partidos, fuera del alcance).\n")

# --- Control de la identidad BS = REL − RES + UNC ---
# La identidad se cumple exactamente solo con bins infinitamente finos.
# Con n_bins = 10, sustituir las probabilidades individuales por el promedio
# del bin introduce una discrepancia pequeña (~0,1%–1,5% del BS). Es esperable
# y no indica error de implementación.
#
# La magnitud de esa discrepancia NO es igual para los tres competidores: es
# mayor cuanto más dispersas son las probabilidades dentro de cada bin. Por
# eso la Sección 5.1 no puede limitarse a restar componentes: tiene que
# reportar qué fracción de la brecha observada explica la descomposición.
cat("\n--- Control de discrepancia por binning ---\n")
descomposicion |>
  transmute(
    fuente,
    discrepancia = round(abs(bs_reconstruido - bs_calculado), 6),
    estado = if_else(abs(bs_reconstruido - bs_calculado) > 0.01, "REVISAR", "OK")
  ) |>
  print(width = Inf)

if (any(abs(descomposicion$bs_reconstruido - descomposicion$bs_calculado) > 0.01)) {
  warning("Discrepancia mayor al 1% del BS en algún competidor — revisar implementación.")
}

# --- Control de UNC ---
# UNC depende solo de los datos observados, así que los tres valores deben
# ser iguales siempre que las tres series se calculen sobre los mismos
# partidos. Si hay faltantes en la serie de cierre, su UNC se calcula sobre
# un subconjunto y puede diferir en los últimos decimales: eso es esperable
# y no indica error.
if (n_distinct(round(descomposicion$incertidumbre, 10)) > 1) {
  if (n_sin_cierre > 0) {
    cat(sprintf(
      "\nUNC difiere levemente en la serie de cierre (%.6f vs. %.6f): esperable,\n",
      descomposicion$incertidumbre[2], descomposicion$incertidumbre[1]))
    cat(sprintf("se calcula sobre %d partidos en lugar de %d.\n",
                nrow(datos_c), nrow(datos)))
  } else {
    warning("UNC difiere entre competidores con cobertura completa: revisar binning.")
  }
} else {
  cat(sprintf("\nUNC idéntica para los tres: %.6f (control OK)\n",
              descomposicion$incertidumbre[1]))
}

cat("\n--- Descomposición por categoría (serie pre-cierre) ---\n")
descomp_mercado$por_categoria |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print()


# -----------------------------------------------------------------------------
# 5.1 DESCOMPOSICIÓN DE LA BRECHA ENTRE COMPETIDORES
# -----------------------------------------------------------------------------
# Las Secciones 4 y 5 responden preguntas distintas y hasta acá quedaban
# desconectadas: la 4 dice SI la diferencia es distinguible de cero, la 5 dice
# cuánto vale cada componente para cada competidor por separado. Falta el
# puente: de la diferencia contrastada en la Sección 4, ¿qué parte viene de
# calibración y qué parte de resolución?
#
# De la identidad BS = REL − RES + UNC, y siendo UNC común a los competidores
# (depende solo de los datos observados), la diferencia entre dos de ellos es:
#
#     ΔBS = ΔREL − ΔRES
#
# Nótese el signo: una PEOR calibración (ΔREL > 0) y una PEOR resolución
# (ΔRES < 0) empujan las dos en la misma dirección, hacia un BS mayor.
#
# POR QUÉ NO ALCANZA CON RESTAR.
# La identidad es exacta solo con bins infinitamente finos. Con 10 bins cada
# competidor arrastra su propia discrepancia de binning, y esas discrepancias
# NO se cancelan al restar. En consecuencia ΔREL − ΔRES no reproduce la
# diferencia observada de Brier Score: la explica solo en parte.
#
# Reportar los componentes sin aclarar esto es un error detectable: cualquiera
# que haga la resta a mano encuentra que los números no cierran. Esta sección
# calcula explícitamente la fracción explicada y deja el residuo a la vista.

descomponer_brecha <- function(nombre_A, nombre_B, fuente_A, fuente_B) {

  comp <- function(f, columna) {
    descomposicion[[columna]][descomposicion$fuente == f]
  }

  d_rel <- comp(fuente_A, "calibracion") - comp(fuente_B, "calibracion")
  d_res <- comp(fuente_A, "resolucion")  - comp(fuente_B, "resolucion")

  # Aportes a ΔBS según la identidad: REL con signo +, RES con signo −
  aporte_rel <-  d_rel
  aporte_res <- -d_res
  explicado  <- aporte_rel + aporte_res

  # Diferencia efectivamente observada y contrastada en la Sección 4
  fila_test <- tests_significacion |>
    filter(metrica == "Brier Score",
           competidor_A == nombre_A, competidor_B == nombre_B)

  observado <- fila_test$dif_media

  # Los porcentajes de reparto solo tienen sentido si ambos aportes empujan
  # en la misma dirección. Si tienen signos opuestos, un componente supera el
  # 100% y el otro es negativo: el número existe pero no se puede leer como
  # "qué parte viene de cada uno". En ese caso se devuelve NA y se explica.
  mismo_signo <- sign(aporte_rel) == sign(aporte_res)

  tibble(
    comparacion      = paste(nombre_A, "-", nombre_B),
    dif_observada    = observado,
    dif_explicada    = explicado,
    pct_explicado    = explicado / observado * 100,
    residuo_binning  = observado - explicado,
    aporte_rel       = aporte_rel,
    aporte_res       = aporte_res,
    pct_rel          = if (mismo_signo) aporte_rel / explicado * 100 else NA_real_,
    pct_res          = if (mismo_signo) aporte_res / explicado * 100 else NA_real_,
    reparto_legible  = mismo_signo,
    significativo    = fila_test$significativo_05,
    p_valor          = fila_test$p_t
  )
}

brecha_descompuesta <- bind_rows(
  descomponer_brecha("Modelo", "Pre-cierre", "Modelo",           "Mercado (pre-cierre)"),
  descomponer_brecha("Modelo", "Cierre",     "Modelo",           "Mercado (cierre)"),
  descomponer_brecha("Cierre", "Pre-cierre", "Mercado (cierre)", "Mercado (pre-cierre)")
)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("5.1 — ¿DE DÓNDE VIENE LA DIFERENCIA? (ΔBS = ΔREL − ΔRES)\n")
cat(strrep("=", 74), "\n", sep = "")

cat("\n--- Cuánto explica la descomposición ---\n")
brecha_descompuesta |>
  select(comparacion, dif_observada, dif_explicada, pct_explicado,
         residuo_binning) |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

cat("\n--- Reparto entre calibración y resolución ---\n")
brecha_descompuesta |>
  select(comparacion, aporte_rel, aporte_res, pct_rel, pct_res,
         significativo) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print(width = Inf)

# --- Lectura por comparación ---
# Se interpreta el reparto SOLO cuando la Sección 4 encontró una diferencia
# distinguible de cero. Descomponer una brecha que no se distingue del azar
# es atribuir origen a un efecto que no está establecido.

cat("\n--- Lectura ---\n")

for (i in seq_len(nrow(brecha_descompuesta))) {
  f <- brecha_descompuesta[i, ]

  cat(sprintf("\n%s (p = %.4f)\n", f$comparacion, f$p_valor))

  if (!f$significativo) {
    cat("  La Sección 4 NO detectó diferencia distinguible de cero.\n")
    cat("  NO interpretar el reparto: se estaría atribuyendo origen a una\n")
    cat("  brecha que no está establecida. Los componentes se reportan solo\n")
    cat("  como descripción, sin lectura causal.\n")
    next
  }

  cat(sprintf("  Diferencia observada: %+.6f | explicada por la descomposición: %+.6f (%.1f%%)\n",
              f$dif_observada, f$dif_explicada, f$pct_explicado))
  cat(sprintf("  Residuo por binning:  %+.6f (no atribuible a REL ni a RES)\n",
              f$residuo_binning))

  if (f$reparto_legible) {
    cat(sprintf("  Calibración (REL): %+.6f  ->  %.1f%% de lo explicado\n",
                f$aporte_rel, f$pct_rel))
    cat(sprintf("  Resolución  (RES): %+.6f  ->  %.1f%% de lo explicado\n",
                f$aporte_res, f$pct_res))

    dominante <- if (abs(f$pct_res) > abs(f$pct_rel)) "RESOLUCIÓN" else "CALIBRACIÓN"
    cat(sprintf("  Componente dominante: %s\n", dominante))

    if (dominante == "RESOLUCIÓN") {
      cat("  REDACTAR: la desventaja no es de honestidad probabilística sino de\n")
      cat("  capacidad de discriminación. El predictor peor no está mal\n")
      cat("  calibrado: se juega menos lejos de la frecuencia base porque no\n")
      cat("  dispone de la información que separa un partido de otro.\n")
    } else {
      cat("  REDACTAR: la desventaja es de calibración. Las probabilidades del\n")
      cat("  predictor peor están sistemáticamente corridas respecto de las\n")
      cat("  frecuencias observadas, con independencia de su discriminación.\n")
    }
  } else {
    cat(sprintf("  Calibración (REL): %+.6f\n", f$aporte_rel))
    cat(sprintf("  Resolución  (RES): %+.6f\n", f$aporte_res))
    cat("  Los dos aportes tienen SIGNOS OPUESTOS: se compensan parcialmente.\n")
    cat("  No corresponde repartir la brecha en porcentajes (uno superaría el\n")
    cat("  100% y el otro sería negativo). Reportar ambos aportes por separado\n")
    cat("  y describir la compensación.\n")
  }
}

cat("\n", strrep("-", 74), "\n", sep = "")
cat("ADVERTENCIA PARA EL TEXTO. Reportá siempre el porcentaje explicado junto\n")
cat("con el reparto. Presentar 'REL aporta X% y RES aporta Y%' sin aclarar\n")
cat("que la descomposición cubre solo una fracción de la brecha invita a que\n")
cat("alguien haga la resta y encuentre que los números no cierran.\n")
cat(strrep("-", 74), "\n", sep = "")


# -----------------------------------------------------------------------------
# 6. MÉTRICAS POR TIPO DE RESULTADO (H, D, A)
# -----------------------------------------------------------------------------
# El Brier Score global mezcla los tres tipos de resultado. Separándolos se
# ve en cuál cada competidor tiene ventaja. Los empates suelen ser los más
# difíciles de predecir para cualquier método.
#
# NOTA DE ESCALA: los bs_* de esta sección son Brier Scores BINARIOS por
# categoría (rango [0,1]), no una desagregación del BS clásico global
# (rango [0,2]). Son un diagnóstico de dónde se concentra el error.
#
# Los ll_* usan log_loss() categórica: al filtrar internamente a ocurrio == 1
# dentro de cada grupo, promedian solo sobre los partidos en que esa
# categoría efectivamente ocurrió (log-loss condicional por tipo de resultado).

metricas_por_resultado <- datos_largo |>
  summarise(
    n          = n(),
    ocurrieron = sum(ocurrio),
    bs_mercado = brier_score(pi_mercado,        ocurrio),
    bs_cierre  = brier_score(pi_mercado_cierre, ocurrio),
    bs_modelo  = brier_score(prob_modelo,       ocurrio),
    ll_mercado = log_loss(pi_mercado,        ocurrio),
    ll_cierre  = log_loss(pi_mercado_cierre, ocurrio),
    ll_modelo  = log_loss(prob_modelo,       ocurrio),
    .by = resultado
  ) |>
  mutate(
    resultado = factor(resultado, levels = c("H", "D", "A"))
  ) |>
  arrange(resultado) |>
  mutate(
    # Mejor competidor en cada fila, por métrica.
    # max.col() sobre la matriz NEGADA devuelve el índice del mínimo.
    mejor_bs = c("Pre-cierre", "Cierre", "Modelo")[
      max.col(-cbind(bs_mercado, bs_cierre, bs_modelo), ties.method = "first")],
    # NO se calcula un "mejor_ll": ver la advertencia de propiedad más abajo.
    # Robustez: ¿cambia algo al pasar de pre-cierre a cierre?
    dif_cierre_vs_precierre = bs_cierre - bs_mercado,
    # Sustantivo: ¿el modelo le gana al mercado en esta categoría?
    dif_modelo_vs_precierre = bs_modelo - bs_mercado
  )

cat("\n", strrep("=", 74), "\n", sep = "")
cat("MÉTRICAS POR TIPO DE RESULTADO — TRES COMPETIDORES\n")
cat(strrep("=", 74), "\n", sep = "")

cat("\n--- Brier Score binario por categoría (menor es mejor) ---\n")
metricas_por_resultado |>
  select(resultado, n, ocurrieron, bs_mercado, bs_cierre, bs_modelo, mejor_bs) |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

cat("\n--- Log-loss condicional por categoría (DESCRIPTIVO, sin ganador) ---\n")
metricas_por_resultado |>
  select(resultado, ocurrieron, ll_mercado, ll_cierre, ll_modelo) |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

cat("\nNO SE DECLARA GANADOR EN ESTA TABLA. log_loss() filtra a ocurrio == 1,\n")
cat("de modo que dentro de una categoría puntúa solo los partidos en que esa\n")
cat("categoría ocurrió. La esperanza de esa cantidad se minimiza reportando\n")
cat("probabilidad 1: premia al predictor que INFLA, esté calibrado o no. No\n")
cat("es una regla de puntuación propia en esta desagregación. Sirve para leer\n")
cat("la MAGNITUD del error condicional, no para ordenar competidores.\n")

cat("\n--- Diferencias (negativo = gana el primer término) ---\n")
metricas_por_resultado |>
  select(resultado, dif_cierre_vs_precierre, dif_modelo_vs_precierre) |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)

cat("\nLECTURA de las dos columnas de diferencias:\n")
cat("  dif_cierre_vs_precierre  -> control de ROBUSTEZ. Si es chica en las\n")
cat("     tres categorías, la conclusión no depende de qué serie se use.\n")
cat("  dif_modelo_vs_precierre  -> resultado SUSTANTIVO. Negativo significa\n")
cat("     que el modelo supera al mercado en esa categoría.\n")
cat("  Ninguna de las dos lleva test: son desagregaciones descriptivas del\n")
cat("     contraste global de la Sección 4. Contrastar las tres por separado\n")
cat("     multiplicaría las comparaciones sin hipótesis pre-especificada.\n")

# Síntesis automática. Solo Brier Score: es la única de las dos métricas que
# conserva la propiedad al desagregar por categoría.
cat("\n--- Síntesis (solo Brier Score) ---\n")
for (i in seq_len(nrow(metricas_por_resultado))) {
  fila <- metricas_por_resultado[i, ]
  cat(sprintf("  %s: mejor en Brier Score = %s\n",
              as.character(fila$resultado), fila$mejor_bs))
}

cat("\nRECORDATORIO (Script 04): el desempeño del modelo en la categoría H debe\n")
cat("leerse junto con el cambio de régimen de la ventaja local declarado en la\n")
cat("Sección 1.1. El modelo tiene la ventaja local congelada en el promedio de\n")
cat("21/22-23/24 y la temporada de validación operó por debajo de ese nivel.\n")


# -----------------------------------------------------------------------------
# 7. MÉTRICAS POR DECIL DE PROBABILIDAD
# -----------------------------------------------------------------------------
# Es la comparación más fina: no solo quién es mejor globalmente, sino quién
# es mejor para los favoritos y quién para los longshots.
#
# Los deciles se definen sobre la probabilidad implícita PRE-CIERRE y esa
# misma partición se usa para los tres competidores. Si cada serie definiera
# sus propios deciles, las columnas no serían comparables entre sí porque
# cada bin agruparía partidos distintos.
#
# -----------------------------------------------------------------------------
# ADVERTENCIA DE INTERPRETACIÓN — EL LOG-LOSS POR DECIL NO ES UNA REGLA PROPIA
#
# La regla general: condicionar sobre el PRONÓSTICO es legítimo; condicionar
# sobre el RESULTADO rompe la propiedad de la regla de puntuación.
#
# Los deciles se definen sobre la probabilidad implícita, es decir sobre el
# pronóstico: eso está bien y el Brier Score por decil es válido, porque se
# calcula sobre las 114 filas del decil ocurran o no.
#
# El log-loss, en cambio, hace las dos cosas a la vez: log_loss() filtra a
# ocurrio == 1, de modo que dentro del decil puntúa SOLO las filas donde el
# evento ocurrió. La esperanza de esa cantidad se minimiza reportando f = 1,
# así que premia al predictor que asigna probabilidades más altas en ese
# decil, esté calibrado o no.
#
# LA DEMOSTRACIÓN ESTÁ EN ESTOS MISMOS DATOS. En el decil 1 la frecuencia
# observada ronda 0,08; el mercado asigna ~0,11 y el modelo ~0,16. El modelo
# está claramente PEOR calibrado ahí — el Brier Score, que es propio, lo dice
# — y sin embargo "gana" el log-loss condicional, precisamente porque infla.
# Misma data, dos respuestas opuestas, y solo una de las dos reglas es válida.
#
# POR ESO: no se calcula mejor_ll y el log-loss no entra en el recuento de
# deciles ganados. La tabla se conserva porque la MAGNITUD del error
# condicional es informativa, pero no sirve para ordenar competidores.
#
# (Comparar log-loss ENTRE deciles tampoco tiene sentido: la submuestra que
# sobrevive al filtro pasa de ~8% de las filas en el decil 1 a ~70% en el 10.
# La columna n_ocurrieron deja visible ese tamaño.)
# -----------------------------------------------------------------------------

metricas_por_bin <- datos_largo |>
  mutate(bin = ntile(pi_mercado, 10)) |>
  summarise(
    n              = n(),
    n_ocurrieron   = sum(ocurrio),
    pi_mercado_med = mean(pi_mercado),
    pi_cierre_med  = mean(pi_mercado_cierre, na.rm = TRUE),
    frec_real      = mean(ocurrio),
    bs_mercado     = brier_score(pi_mercado,        ocurrio),
    bs_cierre      = brier_score(pi_mercado_cierre, ocurrio),
    bs_modelo      = brier_score(prob_modelo,       ocurrio),
    ll_mercado     = log_loss(pi_mercado,        ocurrio),
    ll_cierre      = log_loss(pi_mercado_cierre, ocurrio),
    ll_modelo      = log_loss(prob_modelo,       ocurrio),
    .by = bin
  ) |>
  arrange(bin) |>
  mutate(
    mejor_bs = c("Pre-cierre", "Cierre", "Modelo")[
      max.col(-cbind(bs_mercado, bs_cierre, bs_modelo), ties.method = "first")]
    # Sin mejor_ll: el log-loss por decil no es una regla propia (ver arriba).
  )

cat("\n", strrep("=", 74), "\n", sep = "")
cat("MÉTRICAS POR DECIL — TRES COMPETIDORES\n")
cat(strrep("=", 74), "\n", sep = "")
cat("Deciles definidos sobre la serie pre-cierre (misma partición para los tres)\n")

cat("\n--- Brier Score por decil (menor es mejor) ---\n")
metricas_por_bin |>
  select(bin, n, pi_mercado_med, frec_real,
         bs_mercado, bs_cierre, bs_modelo, mejor_bs) |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(n = 10, width = Inf)

cat("\n--- Log-loss condicional por decil (DESCRIPTIVO, sin ganador) ---\n")
metricas_por_bin |>
  select(bin, n, n_ocurrieron, ll_mercado, ll_cierre, ll_modelo) |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(n = 10, width = Inf)

cat("\nNO se declara ganador: ver la advertencia de propiedad al inicio de la\n")
cat("sección. La columna útil de esta tabla es la MAGNITUD del error\n")
cat("condicional, no el orden entre competidores.\n")

# --- Evidencia numérica de por qué no se declara ganador ---
# Se contrasta, decil por decil, la probabilidad media que asigna cada
# competidor contra la frecuencia observada. Donde un predictor infla y sin
# embargo "gana" el log-loss condicional, queda a la vista el artefacto.
cat("\n--- Probabilidad media asignada vs. frecuencia observada, por decil ---\n")
datos_largo |>
  mutate(bin = ntile(pi_mercado, 10)) |>
  summarise(
    frec_real     = mean(ocurrio),
    media_mercado = mean(pi_mercado),
    media_modelo  = mean(prob_modelo),
    .by = bin
  ) |>
  arrange(bin) |>
  mutate(sesgo_mercado = media_mercado - frec_real,
         sesgo_modelo  = media_modelo  - frec_real,
         across(where(is.numeric), ~ round(., 4))) |>
  print(n = 10, width = Inf)

cat("\nEn los deciles donde el modelo tiene sesgo_modelo mayor que\n")
cat("sesgo_mercado, el modelo está PEOR calibrado; si aun así aparece con\n")
cat("menor log-loss condicional en esa fila, eso es exactamente el artefacto\n")
cat("que motiva no declarar ganador en la tabla anterior.\n")

# --- Recuento de deciles ganados (solo Brier Score) ---
cat("\n--- Deciles ganados por competidor (solo Brier Score) ---\n")
metricas_por_bin |>
  count(competidor = mejor_bs, name = "deciles") |>
  arrange(desc(deciles)) |>
  print(width = Inf)

cat("\nLECTURA. Ganar más deciles no equivale a ganar globalmente: el Brier\n")
cat("Score global pondera por la cantidad de observaciones y por la magnitud\n")
cat("del error, no por cuántos bins se ganaron. Si un competidor gana la\n")
cat("mayoría de los deciles pero pierde en la métrica global, es porque\n")
cat("pierde por mucho en los pocos que cede. Vale la pena comentarlo.\n")

# --- Dónde está la mayor discrepancia entre las dos series de mercado ---
# Control de robustez a nivel decil: si la brecha es chica en todos los
# deciles, la elección de serie no afecta las conclusiones en ningún rango
# de probabilidad, no solo en el promedio.
cat("\n--- Discrepancia entre series de mercado por decil ---\n")
metricas_por_bin |>
  transmute(
    bin,
    pi_mercado_med  = round(pi_mercado_med, 4),
    pi_cierre_med   = round(pi_cierre_med, 4),
    dif_pi          = round(pi_cierre_med - pi_mercado_med, 4),
    dif_bs          = round(bs_cierre - bs_mercado, 6)
  ) |>
  print(n = 10, width = Inf)

cat("\n(dif_bs negativo = la serie de cierre es mejor en ese decil)\n")


# -----------------------------------------------------------------------------
# 8. GRÁFICO 1: BRIER SCORE POR DECIL — TRES COMPETIDORES
# -----------------------------------------------------------------------------

metricas_por_bin |>
  select(bin, bs_mercado, bs_cierre, bs_modelo) |>
  pivot_longer(cols = c(bs_mercado, bs_cierre, bs_modelo),
               names_to = "fuente", values_to = "brier") |>
  mutate(fuente = case_when(
    fuente == "bs_mercado" ~ "Mercado (pre-cierre)",
    fuente == "bs_cierre"  ~ "Mercado (cierre)",
    TRUE                   ~ "Modelo Maher"
  )) |>
  ggplot(aes(x = factor(bin), y = brier, color = fuente, group = fuente)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = etiquetas_series) +
  labs(
    title    = "Brier Score por decil de probabilidad",
    subtitle = "Menor es mejor | Deciles definidos sobre la serie pre-cierre | Premier League 2024/25",
    x        = "Decil de probabilidad implícita (1 = longshots, 10 = favoritos)",
    y        = "Brier Score",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 9. GRÁFICO 2: LOG-LOSS POR DECIL — TRES COMPETIDORES
# -----------------------------------------------------------------------------

metricas_por_bin |>
  select(bin, ll_mercado, ll_cierre, ll_modelo) |>
  pivot_longer(cols = c(ll_mercado, ll_cierre, ll_modelo),
               names_to = "fuente", values_to = "logloss") |>
  mutate(fuente = case_when(
    fuente == "ll_mercado" ~ "Mercado (pre-cierre)",
    fuente == "ll_cierre"  ~ "Mercado (cierre)",
    TRUE                   ~ "Modelo Maher"
  )) |>
  ggplot(aes(x = factor(bin), y = logloss, color = fuente, group = fuente)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = etiquetas_series) +
  labs(
    title    = "Log-loss condicional por decil de probabilidad",
    subtitle = paste0(
      "Premier League 2024/25 | Calculado solo sobre los resultados que ocurrieron\n",
      "DESCRIPTIVO: al condicionar sobre el resultado deja de ser una regla de ",
      "puntuación propia y no ordena competidores"),
    x        = "Decil de probabilidad implícita (1 = longshots, 10 = favoritos)",
    y        = "Log-loss condicional",
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 10. GRÁFICO 3: DESCOMPOSICIÓN DEL BRIER SCORE — TRES COMPETIDORES
# -----------------------------------------------------------------------------
# La incertidumbre (UNC) es igual para los tres y se excluye del gráfico:
# depende solo de los datos, no del predictor.

tibble(
  componente = factor(
    rep(c("Calibración (REL)", "Resolución (RES)"), 3),
    levels = c("Calibración (REL)", "Resolución (RES)")
  ),
  fuente = c(rep("Mercado (pre-cierre)", 2),
             rep("Mercado (cierre)",     2),
             rep("Modelo Maher",         2)),
  valor  = c(
    descomp_mercado$REL, descomp_mercado$RES,
    descomp_cierre$REL,  descomp_cierre$RES,
    descomp_modelo$REL,  descomp_modelo$RES
  )
) |>
  ggplot(aes(x = componente, y = valor, fill = fuente)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_text(aes(label = round(valor, 5)),
            position = position_dodge(width = 0.9),
            vjust = -0.4, size = 2.9) +
  scale_fill_manual(values = etiquetas_series) +
  labs(
    title    = "Descomposición del Brier Score (Murphy, 1973)",
    subtitle = sprintf(
      "REL: menor es mejor (calibración) | RES: mayor es mejor (discriminación)\nIncertidumbre (UNC = %.4f) excluida por ser idéntica para los tres predictores",
      descomp_mercado$UNC
    ),
    x        = NULL,
    y        = "Valor",
    fill     = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 11. GRÁFICO 4: BRIER SCORE POR TIPO DE RESULTADO — TRES COMPETIDORES
# -----------------------------------------------------------------------------

metricas_por_resultado |>
  select(resultado, bs_mercado, bs_cierre, bs_modelo) |>
  pivot_longer(cols = c(bs_mercado, bs_cierre, bs_modelo),
               names_to = "fuente", values_to = "brier") |>
  mutate(
    fuente = case_when(
      fuente == "bs_mercado" ~ "Mercado (pre-cierre)",
      fuente == "bs_cierre"  ~ "Mercado (cierre)",
      TRUE                   ~ "Modelo Maher"
    ),
    resultado = factor(resultado, levels = c("H", "D", "A"),
                       labels = c("Local (H)", "Empate (D)", "Visitante (A)"))
  ) |>
  ggplot(aes(x = resultado, y = brier, fill = fuente)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = etiquetas_series) +
  labs(
    title    = "Brier Score por tipo de resultado",
    subtitle = "Menor es mejor | Premier League 2024/25",
    x        = NULL,
    y        = "Brier Score",
    fill     = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 11.b GRÁFICO 5: DIFERENCIAS CON INTERVALO DE CONFIANZA
# -----------------------------------------------------------------------------
# Traduce la Sección 4 a una figura. Es el gráfico que conviene llevar a la
# defensa: muestra de un vistazo qué diferencias cruzan el cero y cuáles no.
# Se grafica solo el Brier Score porque el log-loss está en otra escala y
# superponerlos haría ilegible el eje.

tests_significacion |>
  filter(metrica == "Brier Score") |>
  mutate(
    # El orden se fija por el par de competidores, no por la cadena de texto:
    # así el gráfico no depende del encoding de la etiqueta.
    orden = case_when(
      competidor_A == "Modelo" & competidor_B == "Pre-cierre" ~ 1,
      competidor_A == "Modelo" & competidor_B == "Cierre"     ~ 2,
      TRUE                                                     ~ 3
    ),
    comparacion = fct_reorder(comparacion, -orden)
  ) |>
  ggplot(aes(x = dif_media, y = comparacion)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_errorbar() con orientation = "y": geom_errorbarh() quedó obsoleta en
  # ggplot2 4.0 y emite un aviso de ciclo de vida en cada corrida.
  geom_errorbar(aes(xmin = ic95_inf, xmax = ic95_sup),
                orientation = "y",
                width = 0.15, linewidth = 0.7, color = "#2c3e50") +
  geom_point(aes(color = significativo_05), size = 3.5) +
  scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#7f8c8d"),
                     labels = c("TRUE" = "p < 0,05", "FALSE" = "p ≥ 0,05")) +
  labs(
    title    = "Diferencias en Brier Score con intervalo de confianza del 95%",
    subtitle = paste0("Contrastes pareados sobre los ", nrow(datos),
                      " partidos de validación 2024/25\n",
                      "Diferencia negativa = el primer competidor predice mejor"),
    x        = "Diferencia media en Brier Score",
    y        = NULL,
    color    = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 12. TABLA RESUMEN FINAL
# -----------------------------------------------------------------------------

cat("\n=== TABLA RESUMEN FINAL ===\n")
cat(sprintf("%-28s %14s %14s %14s\n",
            "", "Merc. pre-cierre", "Merc. cierre", "Modelo"))
cat(sprintf("%-28s %14.6f %14.6f %14.6f\n", "Brier Score global",
            bs_mercado, bs_cierre, bs_modelo))
cat("  Brier Score:  menor es mejor (fórmula clásica, rango [0,2])\n")
cat(sprintf("%-28s %14.4f %14.4f %14.4f\n", "Brier Skill Score",
            bss_mercado, bss_cierre, bss_modelo))
cat("  Skill Score:  mayor es mejor (mejora relativa sobre la climatología)\n")
cat(sprintf("%-28s %14.6f %14.6f %14.6f\n", "Log-loss global",
            ll_mercado, ll_cierre, ll_modelo))
cat("  Log-loss:     menor es mejor (entropía cruzada categórica)\n")
cat(sprintf("%-28s %14.6f %14.6f %14.6f\n", "Calibración BS (REL)",
            descomp_mercado$REL, descomp_cierre$REL, descomp_modelo$REL))
cat("  REL:          menor es mejor (menos sesgo de calibración)\n")
cat(sprintf("%-28s %14.6f %14.6f %14.6f\n", "Resolución BS (RES)",
            descomp_mercado$RES, descomp_cierre$RES, descomp_modelo$RES))
cat("  RES:          mayor es mejor (más capacidad de discriminación)\n")
cat(sprintf("%-28s %14.6f %14.6f %14.6f\n", "Incertidumbre (UNC)",
            descomp_mercado$UNC, descomp_cierre$UNC, descomp_modelo$UNC))
cat("  UNC:          fija para los tres (dificultad intrínseca del problema)\n")
cat(sprintf("\nN: pre-cierre y modelo = %d partidos | cierre = %d partidos\n",
            nrow(datos), nrow(datos_c)))

cat("\n--- Contraste primario (Sección 4) ---\n")
cat(sprintf("Modelo − Mercado pre-cierre en Brier Score: %+.6f\n",
            primario$dif_media))
cat(sprintf("IC 95%% [%+.6f ; %+.6f] | p = %.4f | Wilcoxon p = %.4f\n",
            primario$ic95_inf, primario$ic95_sup,
            primario$p_t, primario$p_wilcoxon))
cat(sprintf("Tamaño de efecto: d = %.3f | el primer competidor predice mejor\n",
            primario$d_cohen))
cat(sprintf("en el %.1f%% de los partidos\n", primario$prop_gana_A * 100))
cat(sprintf("Conclusión: %s\n", primario$mejor))

cat("\n--- Sensibilidad a los ascendidos (Sección 4.6) ---\n")
if (!is.null(primario_sin_asc)) {
  cat(sprintf("Excluyendo %d partidos de ascendidos (n = %d): diferencia %+.6f, p = %.4f\n",
              sum(datos$ascendido), primario_sin_asc$n_pares,
              primario_sin_asc$dif_media, primario_sin_asc$p_t))
  cat(sprintf("Muestra completa (n = %d):                     diferencia %+.6f, p = %.4f\n",
              primario$n_pares, primario$dif_media, primario$p_t))
} else {
  cat("No hay partidos de ascendidos que excluir.\n")
}

cat("\n--- Calibración: contraste contra el piso de ruido (Sección 5) ---\n")
piso_rel |>
  mutate(across(where(is.numeric), ~ round(., 6))) |>
  print(width = Inf)
if (!any(piso_rel$excede_ruido)) {
  cat("Ninguno supera el piso: la brecha de Brier Score no viene de calibración.\n")
}

cat("\n--- Origen de la brecha (Sección 5.1) ---\n")
brecha_primaria <- brecha_descompuesta |>
  filter(comparacion == paste("Modelo", "-", "Pre-cierre"))

cat(sprintf("La descomposición explica el %.1f%% de la diferencia observada.\n",
            brecha_primaria$pct_explicado))
if (isTRUE(brecha_primaria$reparto_legible)) {
  cat(sprintf("De esa parte: %.1f%% resolución | %.1f%% calibración\n",
              brecha_primaria$pct_res, brecha_primaria$pct_rel))
}


# -----------------------------------------------------------------------------
# 13. GUARDAR
# -----------------------------------------------------------------------------
# Se guardan las métricas, los tests de significación, los scores partido a
# partido y los controles de robustez que consume el Script 09.

saveRDS(
  list(
    metricas_globales = tibble(
      fuente      = c("Mercado (pre-cierre)", "Mercado (cierre)", "Modelo"),
      n           = c(nrow(datos), nrow(datos_c), nrow(datos)),
      brier_score = c(bs_mercado,  bs_cierre,  bs_modelo),
      brier_skill = c(bss_mercado, bss_cierre, bss_modelo),
      log_loss    = c(ll_mercado,  ll_cierre,  ll_modelo)
    ),
    descomposicion = tibble(
      fuente        = c("Mercado (pre-cierre)", "Mercado (cierre)", "Modelo"),
      calibracion   = c(descomp_mercado$REL, descomp_cierre$REL, descomp_modelo$REL),
      resolucion    = c(descomp_mercado$RES, descomp_cierre$RES, descomp_modelo$RES),
      incertidumbre = c(descomp_mercado$UNC, descomp_cierre$UNC, descomp_modelo$UNC)
    ),
    metricas_por_resultado = metricas_por_resultado,
    metricas_por_bin       = metricas_por_bin,
    bs_referencia          = bs_ref,
    n_sin_cierre           = n_sin_cierre,

    # --- Nuevo en esta versión ---
    tests_significacion    = tests_significacion,
    contraste_primario     = primario,
    scores_individuales    = scores_individuales,
    brecha_descompuesta    = brecha_descompuesta,

    # --- Controles de robustez ---
    piso_ruido_rel         = piso_rel,
    sesgo_medio_categoria  = sesgo_medio,
    metricas_sin_ascendidos    = metricas_sin_asc,
    contraste_primario_sin_asc = primario_sin_asc,
    brecha_ascendidos          = comparacion_asc,
    equipos_nuevos             = equipos_nuevos
  ),
  "resultados_calibracion.rds"
)

cat("\nArchivo 'resultados_calibracion.rds' guardado.\n")
cat("Incluye ahora 'tests_significacion' (6 contrastes), 'contraste_primario',\n")
cat("'scores_individuales' (380 scores por competidor) y 'brecha_descompuesta'\n")
cat("(origen de cada diferencia: calibración vs. resolución).\n")
cat("Agregados: 'piso_ruido_rel', 'sesgo_medio_categoria',\n")
cat("'metricas_sin_ascendidos', 'contraste_primario_sin_asc' y\n")
cat("'brecha_ascendidos'.\n")

cat("\nPendientes anotados para el Script 09:\n")
cat("  - Usar la brecha de ventaja local (Script 04) para explicar la\n")
cat("    diferencia de RESOLUCIÓN entre modelo y mercado. Los insumos están\n")
cat("    en 'resultados_maher.rds' (ventaja_local_regimen, brecha_ventaja_pct)\n")
cat("    y la evidencia directa en 'sesgo_medio_categoria'.\n")
cat("  - El gráfico g08b del Script 09 reproduce el log-loss por decil: debe\n")
cat("    llevar la misma advertencia de propiedad que la Sección 7 de este\n")
cat("    script, o presentarse sin lectura comparativa.\n")
cat("  - Reportar que la calibración (REL) de los tres competidores es\n")
cat("    indistinguible del piso de ruido: la brecha es toda de resolución.\n")
