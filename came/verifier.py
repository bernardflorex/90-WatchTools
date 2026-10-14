"""Controle de l'installation : versions et presence du noyau compile.

Appele par `lancer.sh --verifier` et `lancer.bat --verifier`. Utile avant de
chercher ailleurs la cause d'un probleme : il dit quelle version de Python
tourne, quelles dependances repondent, et si le noyau compile est actif - sans
lui tout fonctionne, mais une optimisation prend soixante-huit fois plus de
temps.
"""

import importlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("python        :", sys.version.split()[0])
print("plateforme    :", sys.platform)
for nom, attr in (("PySide6", "__version__"),
                  ("numpy", "__version__"),
                  ("numba", "__version__"),
                  ("ezdxf", "__version__")):
    court = nom.split(".")[0]
    try:
        m = importlib.import_module(nom)
        print("%-13s : %s" % (court, getattr(m, attr, "?")))
    except Exception as exc:
        print("%-13s : ABSENT (%s)" % (court, exc.__class__.__name__))

try:
    import came.noyau as N
    print("noyau compile :", "oui" if N.DISPONIBLE else "non (Python pur)")
    from came.chrono import ChronoParams, synthesize_chrono
    r = synthesize_chrono(ChronoParams())
    print("essai de calcul : %s, couple %.5f N.mm, plat %.4f mm"
          % ("ok" if r.ok else "echec", r.couple_coeur, r.plat))
except Exception as exc:
    print("essai de calcul : impossible (%r)" % exc)
