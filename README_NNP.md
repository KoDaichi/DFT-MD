# Calculation Flow

1. 教師データ作成（DFT計算）  
 - atomsk で cif -> xsf  
 - /phase0_2022.01/bin/conv.py で結晶拡大  
 - /phase0 でDFT計算（スキップする場合はtarを解凍）  

2. NNPの作成  
 - /aenet-2.0.3/bin/generate.x で教師データの変換  
 - /aenet-2.0.3/bin/train.x で学習  
 - /aenet-2.0.3/bin/predict.x で予測  

 - /scripts/extract_energy.py (まだしてない)  


Lammpsで予測する場合 (まだしてない)  
