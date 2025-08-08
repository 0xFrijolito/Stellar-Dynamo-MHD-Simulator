@echo off
setlocal enabledelayedexpansion

:: Crear directorio de salida si no existe
if not exist "build" mkdir build

:: Compilador y flags
set CC=gcc
set CFLAGS=-Wextra -g
set INCLUDES=-I./src
set OUTPUT=build/mhd.exe
set OPTIMICE = -O3 

:: Recopilar archivos fuente
set SOURCES=
for /r "src" %%f in (*.c) do (
    set SOURCES=!SOURCES! "%%f"
)

:: Compilar
echo Compilando %OUTPUT%...
%CC% %CFLAGS% %INCLUDES% %SOURCES% -o %OUTPUT% %OPTIMICE%

if %ERRORLEVEL% EQU 0 (
    echo Compilacion exitosa
    echo El ejecutable se encuentra en: %OUTPUT%
) else (
    echo Error durante la compilacion
)

endlocal