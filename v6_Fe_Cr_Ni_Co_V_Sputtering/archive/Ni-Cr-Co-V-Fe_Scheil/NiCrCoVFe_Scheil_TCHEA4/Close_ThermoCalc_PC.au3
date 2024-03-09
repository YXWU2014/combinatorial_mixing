#Include <GuiToolBar.au3>
#include <GuiMenu.au3>
#include <Constants.au3>
;#include "C:\Users\yuxiangw\Documents\Hello_DICTRA\FUNCLib\SysTray_UDF.au3"


; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

; ---0. Kill the Thermo-Calc Program
WinClose("Thermo-Calc 2021a", "")
Sleep(500)

#cs
; ---1. Run the Thermo-Calc Program
Run("C:\Program Files\Thermo-Calc\2016a\Thermo-Calc.exe")
Sleep(1000)

;--------------------------------------------------
; ---2. Activate the Thermo-Calc Window
WinSetState("Thermo-Calc 2016a","",@SW_RESTORE)
DIM $index=_SysTrayIconIndex("Thermo-Calc 2016a", 1)
; MesgBox(1,"index",$index)
_GUICtrlToolbar_ClickIndex(ControlGetHandle('[CLASS:Shell_TrayWnd]','','ToolbarWindow321'), $index, "left",False,2)

WinWaitActive("Thermo-Calc 2016a")
;WinActivate("Thermo-Calc 2016a")
Sleep(100)

; ---3. Resize the Thermo-Calc window to a defined position and size
WinMove("Thermo-Calc 2016a", "", 0, 0, 960, 520)
Sleep(100)

; ---4. Click on the Console window
MouseClick($MOUSE_CLICK_LEFT, 240, 260, 1)
Sleep(100)
#ce

#cs

; ---2. Activate the Thermo-Calc Window
WinSetState("[CLASS:SunAwtFrame]","",@SW_RESTORE)
DIM $index=_SysTrayIconIndex("[CLASS:SunAwtFrame]", 1)
; MesgBox(1,"index",$index)
_GUICtrlToolbar_ClickIndex(ControlGetHandle('[CLASS:Shell_TrayWnd]','','ToolbarWindow321'), $index, "left",False,2)

WinWaitActive("[CLASS:SunAwtFrame]")
; WinActivate("[CLASS:SunAwtFrame]")
Sleep(100)

; ---3. Resize the Thermo-Calc window to a defined position and size
WinMove("[CLASS:SunAwtFrame]", "", 0, 0, 960, 520)
Sleep(100)

; ---4. Click on the Console window
MouseClick($MOUSE_CLICK_LEFT, 240, 260, 1)
Sleep(100)

;Send(@LF & "macro \\ad.monash.edu\home\User072\yuxiangw\Documents\ThermoCalc\terative_Flux-Acr_CemPPT_DISSOLUTION.DCM {ENTER}")

; Sleep(5000)
; WinKill("[CLASS:SunAwtFrame]", "")

#ce

#cs

WinSetState("[CLASS:SunAwtFrame]","",@SW_RESTORE)
DIM $index=_SysTrayIconIndex("[CLASS:SunAwtFrame]", 1)
; MesgBox(1,"index",$index)
_GUICtrlToolbar_ClickIndex(ControlGetHandle('[CLASS:Shell_TrayWnd]','','ToolbarWindow321'), $index, "left",False,2)

; Sleep(3000)
WinMenuSelectItem("[CLASS:SunAwtFrame]", "", "&Window", "Reset &Windows")

; send("come to the hub")
; send("{ENTER}")
; WinSetState("Thermo-Calc 2016a" ,"",@SW_MINIMIZE)


#ce

#comments-start

_Main()

Func _Main()
Run("notepad.exe")
WinWaitActive("[CLASS:Notepad]")
$hWnd = WinGetHandle("[CLASS:Notepad]")
$hMain = _GUICtrlMenu_GetMenu($hWnd)
$t1 = _GUICtrlMenu_GetItemText($hMain, 0)
$hFile = _GUICtrlMenu_GetItemSubMenu($hMain, 0)
$t2 = _GUICtrlMenu_GetItemText($hFile, 3)
Writeln(@LF & "Question: " & @LF & @LF & "How to select submenu entry(""" & $t2 & """) of menu(""" & $t1 & """) ? " & @LF & @LF)

; here is the solution. It is language independent.

WinMenuSelectItem($hWnd, "", $t1, $t2)


EndFunc ;==>_Main

; Write a line of text to Notepad
Func Writeln($sText)
ControlSend("[CLASS:Notepad]", "", "Edit1", $sText & @CR)
EndFunc ;==>Writeln

#comments-end


#comments-start

#Include <GuiToolBar.au3>
#include "SysTray_UDF.au3"
WinSetState("Reliance Netconnect" ,"",@SW_RESTORE)
DIM $index=_SysTrayIconIndex("Reliance Netconnect", 1)
;MsgBox(1,"index",$index)
_GUICtrlToolbar_ClickIndex(ControlGetHandle('[CLASS:Shell_TrayWnd]','','ToolbarWindow321'), $index, "left",False,2)
WinActivate("Reliance Netconnect")
WinWaitActive("Reliance Netconnect")
send("{ENTER}")
send("!c")
WinSetState("Reliance Netconnect" ,"",@SW_MINIMIZE)

#comments-end






