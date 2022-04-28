import PySimpleGUI as sg

layout = [  [sg.Text('My Window')],
            [sg.Input(key='-IN-'), sg.Text(size=(12,1), key='-OUT-')],
            [sg.Button('Go'), sg.Button('Exit')]  ]

window = sg.Window('Window Title', layout, web_port=2222, web_start_browser=False)

while True:             # Event Loop
    event, values = window.read()
    print(event, values)
    if event in (None, 'Exit'):
        break
    if event == 'Go':
        window['-OUT-'].Update(values['-IN-'])
window.close()