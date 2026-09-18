echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v251\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\v251\fluent\ntbin\win64\tell.exe" DESKTOP-PHN5430 63641 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\v251\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="DESKTOP-PHN5430" (%KILL_CMD% 17500) 
if /i "%LOCALHOST%"=="DESKTOP-PHN5430" (%KILL_CMD% 8472) 
if /i "%LOCALHOST%"=="DESKTOP-PHN5430" (%KILL_CMD% 3940)
del "C:\Users\mackm\source\CFDSolver\FLUENT validation files\cleanup-fluent-DESKTOP-PHN5430-8472.bat"
