#!/usr/bin/env python3
"""Traduce una secuencia de nucleótidos a aminoácidos en los 6 marcos de lectura.

Uso:
    python traducir.py ATGGCCATTGTAATGGGCCGC
    python traducir.py -f secuencia.fasta
    echo ATGGCC... | python traducir.py
"""
import argparse
import sys
from itertools import product

BASES = "TCAG"
AMINOACIDOS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_A_AA = {"".join(c): aa for c, aa in zip(product(BASES, repeat=3), AMINOACIDOS)}

COMPLEMENTO = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def limpiar(seq):
    seq = "".join(seq.split()).upper().replace("U", "T")
    invalidos = set(seq) - set("ACGTN")
    if invalidos:
        raise ValueError(f"Caracteres no válidos en la secuencia: {sorted(invalidos)}")
    return seq


def reverso_complementario(seq):
    return seq.translate(COMPLEMENTO)[::-1]


def traducir(seq):
    """Traduce en el marco 0; codones incompletos se ignoran, codones con N dan 'X'."""
    return "".join(
        CODON_A_AA.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3)
    )


def seis_marcos(seq):
    seq = limpiar(seq)
    rc = reverso_complementario(seq)
    resultado = {}
    for i in range(3):
        resultado[f"+{i + 1}"] = traducir(seq[i:])
    for i in range(3):
        resultado[f"-{i + 1}"] = traducir(rc[i:])
    return resultado


def leer_fasta(ruta):
    """Devuelve lista de (nombre, secuencia)."""
    registros, nombre, partes = [], None, []
    with open(ruta) as f:
        for linea in f:
            linea = linea.strip()
            if linea.startswith(">"):
                if nombre is not None:
                    registros.append((nombre, "".join(partes)))
                nombre, partes = linea[1:].strip() or "sin_nombre", []
            elif linea:
                partes.append(linea)
    if nombre is not None:
        registros.append((nombre, "".join(partes)))
    elif partes:
        registros.append(("secuencia", "".join(partes)))
    return registros


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("secuencia", nargs="?", help="secuencia de ADN/ARN (o usa -f / stdin)")
    p.add_argument("-f", "--fasta", help="archivo FASTA de entrada")
    args = p.parse_args()

    if args.fasta:
        registros = leer_fasta(args.fasta)
    elif args.secuencia:
        registros = [("secuencia", args.secuencia)]
    elif not sys.stdin.isatty():
        registros = [("secuencia", sys.stdin.read())]
    else:
        p.error("entrega una secuencia, un archivo FASTA (-f) o usa stdin")

    for nombre, seq in registros:
        print(f">{nombre}")
        try:
            for marco, prot in seis_marcos(seq).items():
                print(f"Marco {marco}: {prot}")
        except ValueError as e:
            sys.exit(f"Error: {e}")
        print()


if __name__ == "__main__":
    main()
