"""Leitura dos mapas de estabilidade ( VTK STRUCTURED_POINTS ASCII, i mais rapido ) e dos campos VTK da SimBoltz."""
import numpy as np


def read_map(path):
    """Devolve a matriz [ j, i ] ( linha = eixo y do mapa ) com 0 = estavel e 1 = instavel."""
    lines = open(path).read().split("\n")
    nx, ny = [int(v) for v in next(l for l in lines if l.startswith("DIMENSIONS")).split()[1:3]]
    k = next(n for n, l in enumerate(lines) if l.startswith("LOOKUP_TABLE"))
    vals = np.array(" ".join(lines[k + 1:]).split(), dtype=float)
    return vals[: nx * ny].reshape(ny, nx)


def read_vectors(path):
    """Campo vetorial ( VECTORS ) de rec_velocity / rec_vorticity: devolve array [ ny, nx, 3 ]."""
    lines = open(path).read().split("\n")
    nx, ny, nz = [int(v) for v in next(l for l in lines if l.startswith("DIMENSIONS")).split()[1:4]]
    k = next(n for n, l in enumerate(lines) if l.startswith("VECTORS"))
    vals = np.array(" ".join(lines[k + 1:]).split(), dtype=float)
    return vals[: 3 * nx * ny * nz].reshape(nz, ny, nx, 3)[0]
