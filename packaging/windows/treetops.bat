@echo off
setlocal EnableExtensions

set "ROOT=%~dp0"
set "PATH=%ROOT%;%PATH%"
set "GDAL_DATA=%ROOT%\share\gdal"
set "PROJ_DATA=%ROOT%\share\proj"
set "PROJ_LIB=%ROOT%\share\proj"

"%ROOT%treetops-cli.exe" %*

endlocal
exit /b %ERRORLEVEL%
