# Install script for the treetops application

#Page license
#Page components
#Page directory
#Page instfiles
#UninstPage uninstConfirm
#UninstPage instfiles

OutFile "dijital_treetops.exe"

# TODO: Put dependencies in a shared dir

InstallDir "$PROGRAMFILES64\dijital.ca\treetops"

Section

	SetOutPath $INSTDIR

	File "C:\boost_1_63_0\lib64-msvc-14.0\boost_date_time-vc140-mt-1_63.dll"
	File "C:\boost_1_63_0\lib64-msvc-14.0\boost_system-vc140-mt-1_63.dll"
	File "C:\boost_1_63_0\lib64-msvc-14.0\boost_filesystem-vc140-mt-1_63.dll"
	File "C:\opt\bin\expat.dll"
	File "C:\opt\bin\freexl.dll"
	File "C:\opt\bin\gdal201.dll"
	File "C:\opt\bin\geos.dll"
	File "C:\opt\bin\geos_c.dll"
	File "C:\opt\bin\iconv.dll"
	File "C:\opt\bin\libcurl.dll"
	File "C:\opt\bin\libeay32.dll"
	File "C:\opt\bin\libmysql.dll"
	File "C:\opt\bin\libpq.dll"
	File "C:\opt\bin\libxml2.dll"
	File "C:\opt\bin\openjp2.dll"
	File "C:\opt\bin\proj.dll"
	File "C:\opt\bin\spatialite.dll"
	File "C:\opt\bin\sqlite3.dll"
	File "C:\opt\bin\ssleay32.dll"
	File "C:\opt\bin\xerces-c_3_1.dll"
	File "C:\opt\bin\zlib1.dll"
	File "C:\Qt\5.7\msvc2015_64\bin\Qt5Core.dll"
	File "C:\Qt\5.7\msvc2015_64\bin\Qt5Widgets.dll"
	File "C:\Qt\5.7\msvc2015_64\bin\Qt5Gui.dll"

	# TODO: Load proj database files/etc. Optionally.

	File "..\makefiles\Release\treetops.dll"
	File "..\makefiles\Release\treetops-app.exe"

	WriteUninstaller "$INSTDIR\uninstall.exe"

	CreateShortcut $SMPROGRAMS\dijital.ca\treetops.lnk $INSTDIR\treetops-app.exe
	CreateShortcut $SMPROGRAMS\dijital.ca\treetops\uninstall.lnk $INSTDIR\uninstall.exe

SectionEnd

Section "Uninstall"

	Delete "$INSTDIR\uninstall.exe"
	Delete "$INSTDIR\*" 
	Delete "$SMPROGRAMS\dijital.ca\treetops.lnk"
	Delete "$SMPROGRAMS\dijital.ca\treetops"

SectionEnd


