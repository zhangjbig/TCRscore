# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))        # 确保能 import 到同目录下的 .pyd/.py

def main(argv=None) -> int:
    import predictTCRscore_core as core         # 编译后会优先加载 _Model_pMHC_core*.pyd
    return core.cli(argv)

if __name__ == "__main__":
    raise SystemExit(main())
