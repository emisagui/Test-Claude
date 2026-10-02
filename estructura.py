#!/usr/bin/env python3
"""Predictor de estructura secundaria (Chou-Fasman simplificado).

Asigna a cada aminoácido: H (hélice α), E (lámina β), T (giro) o C (coil).

Uso:
    python estructura.py MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ
"""
import sys

# Parámetros de Chou-Fasman (x100): Pα, Pβ, Pgiro
PARAMS = {
    "A": (142, 83, 66), "R": (98, 93, 95), "N": (67, 89, 156), "D": (101, 54, 146),
    "C": (70, 119, 119), "Q": (111, 110, 98), "E": (151, 37, 74), "G": (57, 75, 156),
    "H": (100, 87, 95), "I": (108, 160, 47), "L": (121, 130, 59), "K": (114, 74, 101),
    "M": (145, 105, 60), "F": (113, 138, 60), "P": (57, 55, 152), "S": (77, 75, 143),
    "T": (83, 119, 96), "W": (108, 137, 96), "Y": (69, 147, 114), "V": (106, 170, 50),
}
NEUTRO = (100, 100, 100)  # residuos desconocidos (X, etc.)


def _media(valores):
    return sum(valores) / len(valores)


def _regiones(p, ventana, minimo):
    """Busca núcleos (>= `minimo` residuos con P>100 en `ventana`) y los extiende."""
    n = len(p)
    mascara = [False] * n
    for i in range(n - ventana + 1):
        if sum(v > 100 for v in p[i:i + ventana]) >= minimo:
            ini, fin = i, i + ventana  # [ini, fin)
            while ini > 0 and _media(p[ini - 1:ini + 3]) >= 100:
                ini -= 1
            while fin < n and _media(p[fin - 3:fin + 1]) >= 100:
                fin += 1
            for k in range(ini, fin):
                mascara[k] = True
    return mascara


def _tramos(mascara):
    """Tramos contiguos [ini, fin) de True."""
    tramos, ini = [], None
    for i, v in enumerate(mascara + [False]):
        if v and ini is None:
            ini = i
        elif not v and ini is not None:
            tramos.append((ini, i))
            ini = None
    return tramos


def _predecir_segmento(seg):
    n = len(seg)
    if n < 6:
        return "C" * n
    pa, pb, pt = zip(*(PARAMS.get(aa, NEUTRO) for aa in seg))
    hel = _regiones(pa, 6, 4)
    lam = _regiones(pb, 5, 3)

    # Solapamientos: gana la propensión media más alta en el tramo solapado
    for ini, fin in _tramos([h and l for h, l in zip(hel, lam)]):
        if _media(pa[ini:fin]) >= _media(pb[ini:fin]):
            for k in range(ini, fin):
                lam[k] = False
        else:
            for k in range(ini, fin):
                hel[k] = False

    est = ["H" if hel[i] else "E" if lam[i] else "C" for i in range(n)]

    # Giros: tetrapéptidos con Pgiro alta sobre residuos aún sin asignar
    for i in range(n - 3):
        if all(e == "C" for e in est[i:i + 4]):
            t = _media(pt[i:i + 4])
            if t > 100 and t > _media(pa[i:i + 4]) and t > _media(pb[i:i + 4]):
                est[i:i + 4] = ["T"] * 4
    return "".join(est)


def predecir(prot):
    """Devuelve una cadena de igual largo que `prot` con H/E/T/C ('*' se conserva)."""
    prot = prot.upper()
    salida, seg = [], ""
    for aa in prot:
        if aa == "*":
            salida.append(_predecir_segmento(seg) + "*")
            seg = ""
        else:
            seg += aa
    salida.append(_predecir_segmento(seg))
    return "".join(salida)


def main():
    prot = "".join(sys.argv[1:]) or sys.stdin.read()
    prot = "".join(prot.split())
    if not prot:
        sys.exit(__doc__)
    print(prot)
    print(predecir(prot))


if __name__ == "__main__":
    main()
