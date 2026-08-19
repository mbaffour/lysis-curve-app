@echo off
rem ====================================================================
rem  Launch this R/Shiny app locally; it opens in your web browser.
rem  Finds R in this order (first hit wins):
rem    1. R_RSCRIPT env var (point it at your Rscript.exe to force a choice)
rem    2. Rscript on PATH (system R, or an active/base conda with r-base)
rem    3. a normal R install under "C:\Program Files\R\R-*"
rem    4. R inside a conda install (Miniconda3/Anaconda3/miniforge3/mambaforge)
rem  Double-click to run.
rem ====================================================================
setlocal EnableDelayedExpansion
cd /d "%~dp0"
set "RS="
if defined R_RSCRIPT if exist "%R_RSCRIPT%" set "RS=%R_RSCRIPT%"
if not defined RS ( where Rscript >nul 2>nul && set "RS=Rscript" )
if not defined RS (
  for /f "delims=" %%R in ('dir /b /o-n "C:\Program Files\R" 2^>nul ^| findstr /b "R-"') do (
    if not defined RS if exist "C:\Program Files\R\%%R\bin\Rscript.exe" set "RS=C:\Program Files\R\%%R\bin\Rscript.exe"
  )
)
if not defined RS (
  for %%C in (
    "%USERPROFILE%\Miniconda3" "%USERPROFILE%\Anaconda3" "%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge"
    "%ProgramData%\Miniconda3" "%ProgramData%\Anaconda3"
  ) do (
    if not defined RS if exist "%%~C\Scripts\Rscript.exe"          set "RS=%%~C\Scripts\Rscript.exe"
    if not defined RS if exist "%%~C\Library\bin\Rscript.exe"      set "RS=%%~C\Library\bin\Rscript.exe"
    if not defined RS if exist "%%~C\envs\r\Scripts\Rscript.exe"   set "RS=%%~C\envs\r\Scripts\Rscript.exe"
    if not defined RS if exist "%%~C\envs\r\Library\bin\Rscript.exe" set "RS=%%~C\envs\r\Library\bin\Rscript.exe"
  )
)
if not defined RS (
  echo Could not find R ^(Rscript^).
  echo   - Install R from https://cran.r-project.org/  ^(or make sure yours is on PATH^), OR
  echo   - set R_RSCRIPT to your Rscript.exe  ^(e.g. inside a conda env^) and re-run.
  pause
  exit /b 1
)
echo Starting - a browser window will open shortly.
echo (First run installs R packages; this can take a few minutes.)
"%RS%" run_app.R
if errorlevel 1 (
  echo.
  echo The app stopped with an error - see the messages above.
  pause
)
