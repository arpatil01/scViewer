Var CustomLegacyCreated

!macro CustomCodePreInstall
	${If} ${FileExists} "$INSTDIR\App\AppInfo\appinfo.ini"
		ReadINIStr $0 "$INSTDIR\App\AppInfo\appinfo.ini" "Version" "PackageVersion"
		${VersionCompare} $0 "110.0.0.0" $R0
		${If} $R0 == 2
			ReadRegStr $1 HKLM "Software\Microsoft\Windows NT\CurrentVersion" "CurrentBuild"
			${If} $1 < 10000 ;Windows 7/8/8.1
				StrCpy $CustomLegacyCreated true
				${GetParent} $INSTDIR $1
				CreateDirectory "$1\GoogleChromePortableLegacyWin7"
				CopyFiles /SILENT "$INSTDIR\*.*" "$1\GoogleChromePortableLegacyWin7"
				WriteINIStr "$1\GoogleChromePortableLegacyWin7\App\AppInfo\AppInfo.ini" "Details" "AppID" "GoogleChromePortableLegacyWin7"
				WriteINIStr "$1\GoogleChromePortableLegacyWin7\App\AppInfo\AppInfo.ini" "Details" "Name" "Google Chrome Portable (Legacy Win7)"
			${EndIf}
		${EndIf}
	${EndIf}
!macroend

!macro CustomCodePostInstall
	${If} $CustomLegacyCreated == true
		${GetParent} $INSTDIR $1
		CopyFiles /SILENT "$INSTDIR\Other\Source\*.reg" "$1\GoogleChromePortableLegacyWin7\Other\Source"
	${EndIf}
!macroend