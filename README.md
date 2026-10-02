# Test-Claude

Proyectos pequeños y sin dependencias externas.

## Herramientas de secuencias (Python 3)

### `traducir.py` — traducción en los 6 marcos de lectura

Traduce una secuencia de nucleótidos (ADN o ARN) a aminoácidos usando el código genético estándar. Calcula los 3 marcos de la hebra directa (`+1`, `+2`, `+3`) y los 3 del reverso complementario (`-1`, `-2`, `-3`).

- Los codones de parada se muestran como `*`.
- Los codones con `N` se traducen como `X`.
- Los codones incompletos al final se ignoran.
- La `U` se convierte en `T`.
- Si la secuencia tiene caracteres no válidos, termina con un error.

```
python traducir.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
python traducir.py -f secuencias.fasta        # archivo FASTA (admite varias secuencias)
echo ATGGCC... | python traducir.py           # por stdin
python traducir.py -e ATGGCC...               # agrega la estructura secundaria
```

Ejemplo de salida:

```
>secuencia
Marco +1: MAIVMGR*KGAR*
Marco +2: WPL*WAAERVPD
...
```

### `estructura.py` — predictor de estructura secundaria

Predice la estructura secundaria de una secuencia de aminoácidos con el método de Chou-Fasman simplificado. Asigna a cada residuo una letra:

| Letra | Significado |
|-------|-------------|
| `H`   | hélice α    |
| `E`   | lámina β    |
| `T`   | giro        |
| `C`   | coil (sin estructura regular) |

Los `*` (codones de parada) se conservan y dividen la proteína en segmentos que se predicen por separado. Se puede usar solo o desde `traducir.py -e`.

```
python estructura.py MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ
```

**Limitaciones:** Chou-Fasman acierta alrededor de 50–60 % por residuo y tiende a sobrepredecir láminas `E`. Sirve como aproximación didáctica; para resultados fiables usa métodos modernos (PSIPRED, AlphaFold). Además, de los 6 marcos normalmente solo uno corresponde a una proteína real.

## `index.html` — lista de tareas

Lista de tareas en un solo archivo: abre `index.html` en el navegador. Permite agregar, marcar y borrar tareas, que se guardan en `localStorage`.
