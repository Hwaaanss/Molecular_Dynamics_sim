# 입력 파일 배치 안내

```
input/
├── proteins/
│   ├── 6VQN.pdb      # 타겟 단백질 1 (RCSB 에서 다운로드)
│   └── 1VJY.pdb      # 타겟 단백질 2
└── ligand/
    └── drug.sdf      # 약물 1종 (두 시스템이 공유)
```

- 단백질 PDB 는 `./scripts/fetch_inputs.sh` 로 RCSB 에서 자동 다운로드할 수 있다.
- 약물 SDF(`drug.sdf`)는 사용자가 직접 배치한다.
  파일명을 바꾸려면 `LIGAND_SDF` 환경변수로 지정 가능:
  `LIGAND_SDF=input/ligand/myDrug.sdf ./scripts/run_pipeline.sh 6VQN`
