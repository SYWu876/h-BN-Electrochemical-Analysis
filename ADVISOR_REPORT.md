# h-BN Electrochemical Analysis 專案檢查與修正報告

日期：2026-05-05  
本地專案路徑：`D:\h-BN\h-BN-Electrochemical-Analysis`  
目前主要工作分支：`main`，追蹤 `magiccub04/main`  
上游原始專案：`https://github.com/SYWu876/h-BN-Electrochemical-Analysis`  
個人 fork：`https://github.com/magiccub04/h-BN-Electrochemical-Analysis`

## 一、目前總結

本次檢查重點是確認此 GitHub companion archive 是否適合公開發布與回傳 Pull Request。已完成 email metadata、授權、依賴範圍、CI/smoke tests，以及 QAOA 輸出檔名一致性的修正。仍待處理的主要問題是 repository inventory、manifest/checksum 與實際檔案不一致，以及量子分支參數目前仍屬硬編碼結果。

目前 `docs/manifest.csv` 列出 106 個檔案；目前 git 追蹤 90 個檔案。比對後仍有：

- manifest 有列但 repo 沒有的檔案：25 個
- repo 有但 manifest 未列的檔案：9 個

這表示 manifest/checksum 還不能作為可靠的 archive integrity record。

## 二、問題清單與目前狀態

| 編號 | 問題 | 影響 | 目前狀態 |
| --- | --- | --- | --- |
| 1 | README、manifest、checksums 與實際 repo inventory 不一致 | 會影響公開 archive 的可信度，審查者無法依 manifest 驗證檔案完整性 | 尚未修正 |
| 2 | `03_surrogate_qaoa_landscape.py` 原本輸出 `Figure_12*.csv`，但 manifest/既有資料使用 `hBN_*.csv` | README 所說的「重建 EIS outputs」與實際腳本行為不一致 | 已修正，但尚未 commit |
| 3 | `02_quantum_branch_comparison.py` 的 continuous/discrete branch parameters 是硬編碼常數 | 只能重建比較表，不能完整重現 quantum-assisted inference 流程 | 尚未修正，需與老師確認重現性要求 |
| 4 | README 與 CITATION.cff email 不一致 | metadata 與聯絡資訊不一致 | 已修正並 commit |
| 5 | 缺少明確授權 | 公開 GitHub repo 的再利用條件不清楚 | 已修正並 commit |
| 6 | `requirements.txt` 只有 lower bounds | 未來套件大版更新可能造成數值或 CSV 格式差異 | 已修正並 commit |
| 7 | 原本沒有最小測試與 CI | 難以確認腳本是否能在乾淨環境執行 | 已修正並 commit |

補充：先前看到的 `S繚s^alpha` 比較像終端機/編碼顯示問題；目前 repo 內實際內容是 `S·s^alpha`，暫不列為實質資料錯誤。

## 三、已完成的變更

### 已 commit 到 `magiccub04/main`

1. `2a09665 修正email`
   - 將 `CITATION.cff` 的 email 改成 `sywu@gms.ndhu.edu.tw`。
   - README 與 CITATION.cff 的聯絡 email 現已一致。

2. `5a356fc 新增 CI 工作流程、資料授權文件及 MIT 授權，並更新 README 和需求檔案`
   - 新增 `LICENSE`：程式碼採 MIT License。
   - 新增 `DATA_LICENSE.md`：資料、圖表、文件採 CC BY 4.0。
   - 更新 `README.md`：加入 License 區段，將公開發布備註改成 license 已包含、DOI/Zenodo 可另加。
   - 更新 `requirements.txt`：
     - `numpy>=1.24,<3.0`
     - `pandas>=2.0,<4.0`
     - `scipy>=1.10,<2.0`
   - 新增 `requirements-dev.txt`，加入 pytest。
   - 新增 `.github/workflows/ci.yml`，GitHub Actions 會在 Python 3.11/3.12 執行 smoke tests。
   - 新增 `tests/test_package_smoke.py`，測試 email 一致性、script compile，以及三支 EIS scripts 能在暫存專案中產生輸出。

### 目前尚未 commit 的變更

目前工作區還有兩個檔案已修改但尚未 commit：

- `scripts/03_surrogate_qaoa_landscape.py`
- `tests/test_package_smoke.py`

修正內容：

- `03_surrogate_qaoa_landscape.py` 改成輸出正式檔名：
  - `hBN_surrogate_slice.csv`
  - `hBN_qaoa_coarse_landscape.csv`
  - `hBN_qaoa_refined_landscape.csv`
  - `hBN_surrogate_fidelity.csv`
- 移除未文件化的輸出：
  - `Figure_12a_hBN_surrogate_slice.csv`
  - `Figure_12b_hBN_qaoa_coarse_landscape.csv`
  - `Figure_12c_hBN_qaoa_refined_landscape.csv`
  - `Figure_12d_hBN_surrogate_fidelity.csv`
  - `hBN_qaoa_summary.csv`
- smoke test 已同步改成檢查正式檔名存在，並確認舊檔名不再產生。

驗證結果：本地執行 `pytest`，結果為 `3 passed`。

## 四、尚待處理的主要問題

### 1. Manifest 與 checksums 仍需重建

目前 `docs/manifest.csv` 與 `docs/checksums.sha256` 仍是舊 inventory。主要不一致包含：

- manifest 有列但 repo 沒有：
  - `.gitignore`
  - `docs/Note_S6_GitHub_v2.md`
  - `figures/main_text/Figure 1.png` 到 `Figure 8.png`
  - `figures/SI/Figure S1.png` 到 `Figure S13.png`
  - TEM 原始圖的 `.tif/.png` 版本
- repo 有但 manifest 未列：
  - `.github/workflows/ci.yml`
  - `LICENSE`
  - `DATA_LICENSE.md`
  - `requirements-dev.txt`
  - `tests/test_package_smoke.py`
  - TEM 實際存在的 `.jpg` 版本
  - `docs/Note_S9_GitHub_v4.md`

建議下一步要決定：

- 若以目前 repo 實際內容為準：重建 `docs/manifest.csv` 與 `docs/checksums.sha256`。
- 若要保留 README 原本宣稱的完整 figure archive：需補回缺失的 figures、`.tif/.png`、Note_S6 與 `.gitignore`。

### 2. Quantum branch reproducibility 還不完整

`02_quantum_branch_comparison.py` 目前直接寫入 continuous/discrete branch 的參數：

- `p_cont`
- `p_disc`

這可以重建比較表，但不能完整重跑 continuous/discrete quantum-assisted optimization。若老師希望 repo 強調「完整可重現」，就需要補上參數來源、搜尋流程或至少在 README/文件中明確說明這是 archival comparison script。

## 五、建議下一步

1. 先 commit 目前 QAOA 檔名修正：

   ```bash
   git add scripts/03_surrogate_qaoa_landscape.py tests/test_package_smoke.py
   git commit -m "Align QAOA landscape output filenames"
   git push
   ```

2. 接著處理 inventory 問題：

   - 建議以目前 repo 實際內容為準，重建 manifest/checksums。
   - 同步更新 README，避免繼續宣稱包含不存在的 figure exports。

3. 最後與老師確認量子分支要求：

   - 若只作 companion archive，可保留硬編碼參數，但加文件說明來源與限制。
   - 若要完整 reproducibility，需要新增 continuous/discrete inference 的重跑流程。

## 六、可向老師簡短報告的版本

目前我已經把專案中比較明顯的發布品質問題先修掉，包括 email metadata 一致、加入 MIT 與 CC BY 4.0 授權、限制 Python 套件版本範圍，並加入 GitHub Actions 與 smoke tests。後續也修了 QAOA landscape 腳本的輸出檔名，讓它與 manifest/既有資料檔名一致。

目前仍待處理的是 archive inventory 問題：manifest/checksum 仍列出不存在的 figures、TEM `.tif/.png` 和舊版 Note_S6，而 repo 實際有 `.jpg` 與 Note_S9。另一個需要討論的是 quantum branch 的可重現性，目前 continuous/discrete branch 參數仍是硬編碼，適合作為比較表重建，但還不是完整 quantum-assisted inference 重跑流程。
