# Tarificación bajo incertidumbre en apuestas deportivas

Código, datos y texto de la tesina *Tarificación bajo incertidumbre en apuestas deportivas: margen, sesgos y calibración frente a un modelo actuarial de goles*, de Agustín Valente, presentada en Práctica Profesional del Actuario (Cód. 728), Facultad de Ciencias Económicas, Universidad de Buenos Aires. Tutor: Gonzalo Ramos. Defendida el 19 de octubre de 2026.

📄 **[Leer la tesina (PDF)](tesina_valente.pdf)** · 📊 **[Ver la presentación de la defensa (PDF)](presentacion_defensa.pdf)**

## De qué trata

Una cuota deportiva se construye como una prima de seguro: una estimación de probabilidad (la prima pura) más un recargo que le asegura a la casa de apuestas un resultado esperado positivo (el *overround*). El trabajo audita esa tarificación sobre la Premier League inglesa y las cuotas de Bet365, en tres pasos:

1. Estima un modelo de Poisson de goles (Maher, 1982) y contrasta dos extensiones: la corrección de Dixon y Coles (1997) para marcadores bajos y un rating Elo como covariable dinámica (Hvattum y Arntzen, 2010).
2. Mide el margen comercial de las cuotas y contrasta la presencia de *favourite-longshot bias*.
3. Compara la calidad de las probabilidades del modelo con la de las implícitas en las cuotas mediante Brier Score, log-loss y la descomposición de Murphy (1973).

## Resultados principales

- **Margen:** el *overround* promedia 5,43 % sobre los 1520 encuentros y se mantiene estable entre temporadas, como una política de precios más que como un recargo por incertidumbre.
- **Modelo:** ninguna de las dos extensiones se retiene. La corrección de Dixon y Coles no aporta ajuste significativo en esta muestra y el rating Elo resulta redundante con los parámetros de ataque y defensa, de modo que la especificación final es la de Maher.
- **Sesgo:** los tres instrumentos aplicados apuntan en la dirección del *favourite-longshot bias*, pero ninguno lo distingue del azar; el diseño sí permite descartar un sesgo de magnitud sustancial.
- **Modelo contra mercado:** el mercado predice mejor sobre la temporada 2024/25 y la diferencia no es atribuible al azar, pero ambos resultan igualmente bien calibrados. La brecha está en la resolución: el mercado discrimina mejor entre partidos porque incorpora información (lesiones, sanciones, cambios de entrenador) que un modelo basado solo en goles no observa.

El detalle de cada resultado, con sus contrastes y salvedades, está en los capítulos 4 y 5 de la tesina.

## Contenido del repositorio

```
01_carga_datos.R … 09_sintesis_final.R   nueve programas en R, en orden de ejecución
E0_2122.csv … E0_2425.csv                 datos de las cuatro temporadas
tesina_valente.pdf                        texto completo de la tesina
presentacion_defensa.pdf                  diapositivas de la defensa oral
CITATION.cff                              cómo citar este trabajo
LICENSE                                   licencias del código, el texto y los datos
```

## Datos

Cuatro temporadas consecutivas de la Premier League, 1520 encuentros en total, publicadas por [Football-Data.co.uk](https://www.football-data.co.uk/englandm.php) y recuperadas en abril de 2026. Los archivos se incluyen sin modificaciones respecto del original.

| Archivo | Temporada | Encuentros | Uso |
|---|---|---|---|
| `E0_2122.csv` | 2021/22 | 380 | Entrenamiento |
| `E0_2223.csv` | 2022/23 | 380 | Entrenamiento |
| `E0_2324.csv` | 2023/24 | 380 | Entrenamiento |
| `E0_2425.csv` | 2024/25 | 380 | Validación |

De las varias decenas de columnas que publica el repositorio se emplean doce: fecha, equipos, goles al final del encuentro, resultado, y las cuotas decimales de Bet365 para el mercado 1X2 en sus dos series, la relevada varios días antes del encuentro (`B365H/D/A`) y la de cierre (`B365CH/CD/CA`).

El análisis del mercado (margen y sesgo) usa las cuatro temporadas, porque no involucra al modelo. La comparación entre el modelo y el mercado usa solo 2024/25, porque debe hacerse fuera de la muestra de estimación.

## Programas

Los nueve programas se ejecutan **en orden** y se comunican mediante archivos intermedios `.rds`: cada uno lee los objetos que produjo el anterior y guarda los propios.

| Programa | Qué hace | Tiempo aproximado |
|---|---|---|
| `01_carga_datos.R` | Carga las cuatro temporadas, valida la integridad de los datos y define la partición entre entrenamiento y validación | segundos |
| `02_poisson_biv_con_correccion.R` | Estima por máxima verosimilitud el modelo de Maher con la corrección τ de Dixon y Coles | 1 a 3 minutos |
| `03_test_rho.R` | Contrasta la significatividad del parámetro de la corrección y construye su verosimilitud perfilada | 20 a 25 minutos |
| `04_modelo_maher.R` | Consolida el modelo adoptado, resuelve la imputación de los equipos ascendidos y genera las probabilidades 1X2 de la temporada de validación | menos de un minuto |
| `05_ratings_elo.R` | Construye el rating Elo, calibra sus parámetros y descompone el diferencial en sus componentes estática y dinámica | segundos |
| `06_modelo_maher_elo.R` | Compara cuatro especificaciones con y sin el diferencial de Elo y selecciona la definitiva por validación temporal interna | segundos |
| `07_overround_flb.R` | Calcula el *overround* y contrasta la presencia de *favourite-longshot bias* | 1 a 2 minutos |
| `08_calibracion.R` | Compara la calidad predictiva del modelo con la del mercado y descompone la diferencia | menos de un minuto |
| `09_sintesis_final.R` | Consolida los resultados en tablas `.csv` y figuras `.png` | menos de un minuto |

## Cómo ejecutar

1. Descargar el repositorio completo (botón verde **Code → Download ZIP**, o `git clone`).
2. Instalar los paquetes necesarios, una sola vez:
   ```r
   install.packages(c("tidyverse", "janitor", "ggrepel"))
   ```
3. Abrir R o RStudio con la carpeta del repositorio como directorio de trabajo. Los programas leen y escriben con rutas relativas, de modo que todos los archivos intermedios se generan allí.
4. Ejecutar los nueve programas en orden numérico, por ejemplo:
   ```r
   for (f in sort(list.files(pattern = "^0[1-9]_.*\\.R$"))) source(f)
   ```

Los archivos `.rds` intermedios, las tablas y las figuras no se incluyen en el repositorio: se generan al ejecutar los programas. El análisis completo tarda alrededor de media hora, dominada por el programa 03.

## Entorno

- R 4.5.3 (2026-03-11 ucrt)
- `tidyverse` 2.0.0
- `janitor`, empleado únicamente en el programa 01 (`clean_names()`)
- `ggrepel`, empleado únicamente en el programa 09 (etiquetas del gráfico de parámetros)

No se utilizan paquetes econométricos externos: la estimación por máxima verosimilitud se resuelve con el optimizador de base del lenguaje y los estimadores de varianza agrupada están implementados de forma explícita.

Todas las simulaciones parten de la semilla 904158, de modo que sus resultados son exactamente reproducibles.

## Verificaciones

Cada programa incorpora controles con umbrales explícitos que detienen la ejecución cuando fallan: de integridad de los datos, de coherencia entre magnitudes que deben coincidir por construcción, y de coherencia con valores conocidos de antemano, como la cantidad de partidos por temporada. Si un programa se interrumpe con un error, el mensaje indica cuál de esas condiciones no se cumplió.

Cada programa puede ejecutarse también en una sesión nueva de R, sin depender de objetos que hayan quedado en memoria de una ejecución anterior. El análisis completo se reprodujo de esa forma desde cero en octubre de 2026 (R 4.3.3, `tidyverse` 2.0.0) y arroja las mismas cifras que reporta la tesina.

## Cómo citar

> Valente, A. (2026). *Tarificación bajo incertidumbre en apuestas deportivas: margen, sesgos y calibración frente a un modelo actuarial de goles* [Tesina de grado, Facultad de Ciencias Económicas, Universidad de Buenos Aires]. https://github.com/agusvalente/tesina-tarificacion-apuestas-deportivas

GitHub también ofrece la cita en otros formatos desde el botón **Cite this repository** de la barra lateral.

## Licencia

- Código (archivos `.R`): [MIT](LICENSE).
- Texto de la tesina y presentación (`tesina_valente.pdf`, `presentacion_defensa.pdf`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.es). Se puede compartir y adaptar citando la autoría. Los logos institucionales que figuran en la presentación pertenecen a la Universidad de Buenos Aires y no están cubiertos por esta licencia.
- Datos (`E0_*.csv`): pertenecen a [Football-Data.co.uk](https://www.football-data.co.uk/) y no están cubiertos por las licencias anteriores.
