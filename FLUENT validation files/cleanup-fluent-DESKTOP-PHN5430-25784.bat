echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v251\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\v251\fluent\ntbin\win64\tell.exe" DESKTOP-PHN5430 56641 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\v251\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="DESKTOP-PHN5430" (%KILL_CMD% 2076) 
if /i "%LOCALHOST%"=="DESKTOP-PHN5430" (%KILL_CMD% 25784) 
if /i "%LOCALHOST%"=="DESKTOP-PHN5430" (%KILL_CMD% 29144)
del "C:\Users\mackm\source\CFDSolver\FLUENT validation files\cleanup-fluent-DESKTOP-PHN5430-25784.bat"
