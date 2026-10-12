"""Icones des outils, dessinees au trait dans les couleurs du theme.

Une icone est une fonction `f(q, c)` qui dessine dans un carre de cote 1,
centre sur l'origine, y vers le haut ; `c` est sa couleur.
"""

import math

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import (QBrush, QColor, QIcon, QPainter, QPainterPath, QPen,
                           QPixmap, QPolygonF)


def _polygone(points):
    return QPolygonF([QPointF(x, y) for x, y in points])


def _aplat(c):
    f = QColor(c)
    f.setAlpha(70)
    return QBrush(f)


def came(q, c):
    """Coeur de remise a zero : le rayon croit lineairement avec l'angle,
    du creux (en haut, pres de l'axe) a la pointe (en bas)."""
    n, r0, r1, cy = 120, 0.10, 0.56, 0.20
    pts = []
    for i in range(n + 1):
        t = -math.pi + 2.0 * math.pi * i / n          # 0 : le creux
        r = r0 + (r1 - r0) * abs(t) / math.pi
        pts.append((r * math.sin(t), cy + r * math.cos(t)))
    q.setBrush(_aplat(c))
    q.drawPolygon(_polygone(pts))
    q.setBrush(QBrush(c))
    q.drawEllipse(QPointF(0.0, cy), 0.035, 0.035)


def engrenage(q, c):
    """Roue dentee."""
    z, pied, tete = 12, 0.31, 0.43
    pas = 2.0 * math.pi / z
    pts = []
    for k in range(z):
        a = k * pas
        for da, r in ((-0.27, pied), (-0.15, tete), (0.15, tete), (0.27, pied)):
            pts.append((r * math.cos(a + da * pas * 2), r * math.sin(a + da * pas * 2)))
    roue = QPainterPath()
    roue.addPolygon(_polygone(pts))
    roue.closeSubpath()
    roue.addEllipse(QPointF(0, 0), 0.11, 0.11)       # alesage : reste vide
    q.setBrush(_aplat(c))
    q.drawPath(roue)


def spiral(q, c):
    """Spiral d'Archimede, de la virole au piton."""
    tours, r0, r1, n = 4.0, 0.07, 0.43, 400
    pts = []
    for i in range(n + 1):
        a = 2.0 * math.pi * tours * i / n
        r = r0 + (r1 - r0) * i / n
        pts.append((r * math.cos(a), r * math.sin(a)))
    q.setBrush(Qt.BrushStyle.NoBrush)
    q.drawPolyline(_polygone(pts))
    q.setBrush(QBrush(c))
    q.drawEllipse(QPointF(0, 0), 0.03, 0.03)
    q.drawEllipse(QPointF(*pts[-1]), 0.03, 0.03)


def icone(dessin, couleur, cote=128):
    """QIcon de `cote` pixels a partir d'une fonction de dessin."""
    pm = QPixmap(cote, cote)
    pm.fill(Qt.GlobalColor.transparent)
    q = QPainter(pm)
    q.setRenderHint(QPainter.RenderHint.Antialiasing)
    q.translate(cote / 2.0, cote / 2.0)
    q.scale(cote, -cote)
    crayon = QPen(couleur, 0.035)
    crayon.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    crayon.setCapStyle(Qt.PenCapStyle.RoundCap)
    q.setPen(crayon)
    dessin(q, couleur)
    q.end()
    return QIcon(pm)
