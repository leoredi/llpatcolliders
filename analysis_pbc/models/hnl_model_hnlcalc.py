
from __future__ import annotations

import ast
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, Any


_ALLOWED_HNLCALC_EXPR_CALLS = {
    "get_2body_br",
    "get_2body_br_tau",
    "get_3body_dbr_pseudoscalar",
    "get_3body_dbr_vector",
    "get_3body_dbr_baryon",
    "get_3body_dbr_tau",
}

_ALLOWED_NUMPY_ATTRS = {
    "sqrt",
    "pi",
}


@lru_cache(maxsize=4096)
def _parse_safe_expr(expr: str) -> ast.AST:
    return ast.parse(expr, mode="eval").body


class _SafeExprEvaluator:
    def __init__(self, *, hnl: Any, mass: float, coupling: float, np_module: Any) -> None:
        self._hnl = hnl
        self._mass = float(mass)
        self._coupling = float(coupling)
        self._np = np_module

    def eval(self, expr: str):
        return self._eval_node(_parse_safe_expr(expr))

    def _eval_node(self, node: ast.AST):
        if isinstance(node, ast.Constant):
            if isinstance(node.value, (int, float, str)):
                return node.value
            raise ValueError(f"Unsupported constant type: {type(node.value).__name__}")

        if isinstance(node, ast.Name):
            if node.id == "hnl":
                return self._hnl
            if node.id == "mass":
                return self._mass
            if node.id == "coupling":
                return self._coupling
            if node.id == "np":
                return self._np
            raise ValueError(f"Name not allowed in expression: {node.id}")

        if isinstance(node, ast.UnaryOp):
            operand = self._eval_node(node.operand)
            if isinstance(node.op, ast.UAdd):
                return +operand
            if isinstance(node.op, ast.USub):
                return -operand
            raise ValueError(f"Unary operator not allowed: {type(node.op).__name__}")

        if isinstance(node, ast.BinOp):
            left = self._eval_node(node.left)
            right = self._eval_node(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                return left**right
            raise ValueError(f"Binary operator not allowed: {type(node.op).__name__}")

        if isinstance(node, ast.Attribute):
            if not isinstance(node.value, ast.Name):
                raise ValueError("Only simple attribute access is allowed (e.g. hnl.foo, np.sqrt)")
            base = node.value.id
            attr = node.attr
            if attr.startswith("_"):
                raise ValueError("Access to private/dunder attributes is not allowed")
            if base == "hnl":
                if attr not in _ALLOWED_HNLCALC_EXPR_CALLS:
                    raise ValueError(f"hnl.{attr} is not allowed in expressions")
                return getattr(self._hnl, attr)
            if base == "np":
                if attr not in _ALLOWED_NUMPY_ATTRS:
                    raise ValueError(f"np.{attr} is not allowed in expressions")
                return getattr(self._np, attr)
            raise ValueError(f"Attribute access not allowed on base: {base}")

        if isinstance(node, ast.Call):
            if node.keywords:
                raise ValueError("Keyword arguments are not allowed in expressions")
            func = self._eval_node(node.func)
            args = [self._eval_node(a) for a in node.args]
            return func(*args)

        raise ValueError(f"Expression node not allowed: {type(node).__name__}")


THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[1]
HNL_CALC_DIR = REPO_ROOT / "HNLCalc"

if str(HNL_CALC_DIR) not in sys.path:
    sys.path.insert(0, str(HNL_CALC_DIR))

try:
    import HNLCalc as _hnlcalc_module
    if hasattr(_hnlcalc_module, "MODEL_DIR"):
        _hnlcalc_module.MODEL_DIR = HNL_CALC_DIR.parent / "model"
except ImportError as exc:
    raise ImportError(
        "Could not import the local HNLCalc package. "
        "Make sure analysis_pbc/HNLCalc exists and is a valid Python module."
    ) from exc


class HNLModel:

    def __init__(
        self,
        mass_GeV: float,
        Ue2: float,
        Umu2: float,
        Utau2: float,
        **extra_args: Any,
    ) -> None:
        self.mass_GeV = float(mass_GeV)
        self.Ue2 = float(Ue2)
        self.Umu2 = float(Umu2)
        self.Utau2 = float(Utau2)
        self.extra_args = extra_args

        self._hnlcalc = self._build_hnlcalc()


    def _build_hnlcalc(self):
        from HNLCalc import HNLCalc as _HNLCalcClass
        import numpy as np

        ve = np.sqrt(self.Ue2)
        vmu = np.sqrt(self.Umu2)
        vtau = np.sqrt(self.Utau2)

        hnl = _HNLCalcClass(ve=ve, vmu=vmu, vtau=vtau)

        epsilon = np.sqrt(self.Ue2 + self.Umu2 + self.Utau2)

        hnl.get_br_and_ctau(mpts=np.array([self.mass_GeV]), coupling=epsilon)

        return hnl


    @property
    def ctau0_m(self) -> float:
        return self._hnlcalc.ctau[0]


    def production_brs(self) -> Dict[int, float]:
        import numpy as np

        channels_2body = sum(self._hnlcalc.get_channels_2body()["mode"].values(), [])
        channels_3body = sum(self._hnlcalc.get_channels_3body()["mode"].values(), [])

        mass = self.mass_GeV
        coupling = np.sqrt(self.Ue2 + self.Umu2 + self.Utau2)
        hnl = self._hnlcalc

        evaluator = _SafeExprEvaluator(hnl=hnl, mass=mass, coupling=coupling, np_module=np)

        def _safe_eval(expr: str):
            return evaluator.eval(expr)

        br_per_parent = {}

        for channel in channels_2body:
            pid0 = int(channel['pid0'])
            pid1 = int(channel['pid1'])
            br_string = channel['br']

            m_parent = hnl.masses(str(pid0))
            m_daughter = hnl.masses(str(pid1))

            if mass < m_parent - m_daughter:
                br_formula = _safe_eval(br_string)
                br_value = _safe_eval(br_formula)

                parent_abs = abs(pid0)
                br_per_parent[parent_abs] = br_per_parent.get(parent_abs, 0.0) + br_value

        for channel in channels_3body:
            pid0 = int(channel['pid0'])
            pid1 = int(channel['pid1'])
            pid2 = int(channel['pid2'])
            br_string = channel['br']
            integration = channel['integration']

            m0 = hnl.masses(str(pid0))
            m1 = hnl.masses(str(pid1))
            m2 = hnl.masses(str(pid2))

            if mass < m0 - m1 - m2:
                br_diff = _safe_eval(br_string)

                br_value = hnl.integrate_3body_br(
                    br_diff, mass, m0, m1, m2,
                    coupling=coupling,
                    integration=integration
                )

                parent_abs = abs(pid0)
                br_per_parent[parent_abs] = br_per_parent.get(parent_abs, 0.0) + br_value


        m_W = 80.4
        m_Z = 91.2

        if mass < m_W:
            BR_W_to_lnu_SM = 0.1086

            r_W = mass / m_W
            phase_space_W = (1.0 - r_W**2)**2
            helicity_W = (1.0 + r_W**2)

            br_W = (self.Ue2 + self.Umu2 + self.Utau2) * BR_W_to_lnu_SM * phase_space_W * helicity_W
            br_per_parent[24] = br_W

        if mass < m_Z:
            BR_Z_to_nunu_SM = 0.201 / 3.0

            r_Z = mass / m_Z
            phase_space_Z = (1.0 - r_Z**2)**2
            helicity_Z = (1.0 + r_Z**2)

            br_Z = (self.Ue2 + self.Umu2 + self.Utau2) * BR_Z_to_nunu_SM * phase_space_Z * helicity_Z
            br_per_parent[23] = br_Z

        return br_per_parent


    def __repr__(self) -> str:
        return (
            f"HNLModel(mass_GeV={self.mass_GeV}, "
            f"Ue2={self.Ue2}, Umu2={self.Umu2}, Utau2={self.Utau2})"
        )
