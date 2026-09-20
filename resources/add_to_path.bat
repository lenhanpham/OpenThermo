@echo off
REM Script to add current directory to user PATH (no admin required)

echo Adding Rotbond to User PATH...
echo Current directory: %~dp0
echo.

REM Get current directory without trailing backslash
set "CURRENT_DIR=%~dp0"
set "CURRENT_DIR=%CURRENT_DIR:~0,-1%"

REM Check if directory is already in user PATH
for /f "usebackq tokens=2*" %%A in (`reg query "HKEY_CURRENT_USER\Environment" /v PATH 2^>nul`) do set "USER_PATH=%%B"

REM Handle case where user PATH doesn't exist
if "%USER_PATH%"=="" (
    set "NEW_PATH=%CURRENT_DIR%"
    echo Creating new user PATH variable...
) else (
    REM Check if already in PATH
    echo %USER_PATH% | find /i "%CURRENT_DIR%" >nul
    if not errorlevel 1 (
        echo Directory is already in user PATH!
        echo No changes needed.
        goto :end
    )
    set "NEW_PATH=%USER_PATH%;%CURRENT_DIR%"
)

REM Update user PATH in registry
reg add "HKEY_CURRENT_USER\Environment" /v PATH /t REG_EXPAND_SZ /d "%NEW_PATH%" /f

if errorlevel 1 (
    echo Error: Failed to update user PATH
    pause
    exit /b 1
)

REM Broadcast environment change to current session
setx PATH "%NEW_PATH%" >nul

echo.
echo SUCCESS: Added to user PATH!
echo Path added: %CURRENT_DIR%
echo.
echo The change will take effect for new command prompts and applications.
echo Current command prompt may need to be restarted to see the change.

:end
echo.
echo You can now run 'OpenThermo' from any directory!
echo.
pause