from tkinter import *
from tkinter import ttk
from tkinter import messagebox
root=Tk()
root.title('Физический калькулятор')
root.geometry('500x300')
Melt_T = 0
clicks = 0
WindowCount = 0
S_Mass = []
Sum_Hc = 0
Energy = 0
Melt_Data = []
Melt_Data_H = []
Temperatur = []



tab_control = ttk.Notebook(root)

Menu = ttk.Frame(tab_control)
tab_control.add(Menu, text='Меню')

def Btn_take1():
    global clicks
    global Mass
    global H_capacity
    global Temperature
    global S_Mass
    global AverageHc
    global Energy
    global Sum_Hc
    global Melt_T
    global Melt_H
    global Melt_Data
    global Melt_Data_H
    global Temperatur
    S_Mass.append(float(Mass.get()))
    Sum_Hc+=float(H_capacity.get())*float(Mass.get())
    Energy+=float(Temperature.get())*float(H_capacity.get())*float(Mass.get())
    clicks+=1
    Melt_Data_H.append(float(Melt_H.get()))
    Melt_Data.append(float(Melt_T.get()))
    Temperatur.append(float(Temperature.get()))
    


def Btn_Thermal():
    global WindowCount
    global clicks
    global Mass
    global H_capacity
    global Temperature
    global Melt_T
    global Melt_H
    WindowCount+=1
    if WindowCount > 1: 
        tab_control.forget(1)
        WindowCount = 1
    Thermal = ttk.Frame(tab_control)
    tab_control.add(Thermal, text='Задачи на тепло')
    Btn_take = Button(Thermal, text="Смешать", command = Btn_take1)
    Btn_take.grid(row=0,column=1)
    message = StringVar()
    Mass = Entry(Thermal)
    Mass.grid(row=1,column=2, padx=5, pady=5)
    lb_m = Label(Thermal, text='Масса(m)')
    lb_m.grid(row=1,column=1)
    H_capacity = Entry(Thermal)
    H_capacity.grid(row=2,column=2, padx=5, pady=5)
    lb_c = Label(Thermal, text='Теплоемкость(с)')
    lb_c.grid(row=2,column=1)
    Temperature = Entry(Thermal)
    Temperature.grid(row=3,column=2, padx=5, pady=5)
    lb_t = Label(Thermal, text='температура(t)')
    lb_t.grid(row=3,column=1)
    Btn_A = Button(Thermal, text='Очистка', command = Clear)
    Btn_A.grid(column=0, row=2)
    Btn_A = Button(Thermal, text='Расчет', command = Calc)
    Btn_A.grid(column=0, row=1)
    Melt_T = Entry(Thermal)
    Melt_T.grid(row=4,column=2, padx=5, pady=5)
    lb_Mt = Label(Thermal, text='Температура плавления')
    lb_Mt.grid(row=4,column=1)
    Melt_H = Entry(Thermal)
    Melt_H.grid(row=5,column=2, padx=5, pady=5)
    lb_Mh = Label(Thermal, text='Теплота плавления')
    lb_Mh.grid(row=5,column=1)


def Btn_С():
    global WindowCount
    global Temperature
    global Temperature0
    global Temperature1
    global Force
    global Power
    global Speed
    global Eff
    global Amperage
    global Voltage
    global Mass
    global CapcTemp
    global Mt
    global Time
    global CapcTemp1
    global r
    WindowCount+=1
    if WindowCount > 1: 
        tab_control.forget(1)
        WindowCount = 1
    Thermal = ttk.Frame(tab_control)
    tab_control.add(Thermal, text='Калькулятор')
    Temperature0 = Entry(Thermal)
    Temperature0.grid(row=1,column=0, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Температура в начале(t0 °С)')
    lb_t0.grid(row=0,column=0)
    Temperature1 = Entry(Thermal)
    Temperature1.grid(row=3,column=0, padx=5, pady=5)
    lb_t1 = Label(Thermal, text='Температура в конце(t1 °С)')
    lb_t1.grid(row=2,column=0)
    Temperature = Entry(Thermal)
    Temperature.grid(row=8,column=0, padx=5, pady=5)
    lb_t1 = Label(Thermal, text='Изменение температуры(Δt °С)')
    lb_t1.grid(row=7,column=0)
    lb_t1 = Label(Thermal, text='Или')
    lb_t1.grid(row=5,column=0)
    Force = Entry(Thermal)
    Force.grid(row=1,column=1, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Сила тяги(F H)')
    lb_t0.grid(row=0,column=1)
    Speed = Entry(Thermal)
    Speed.grid(row=3,column=1, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Скорость(V м/с)')
    lb_t0.grid(row=2,column=1)
    Voltage = Entry(Thermal)
    Voltage.grid(row=5,column=1, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Напряжение(U В)')
    lb_t0.grid(row=4,column=1)
    Amperage = Entry(Thermal)
    Amperage.grid(row=7,column=1, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Сила тока(I А)')
    lb_t0.grid(row=6,column=1)
    Eff = Entry(Thermal)
    Eff.grid(row=1,column=2, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='КПД(η %)')
    lb_t0.grid(row=0,column=2)
    Time = Entry(Thermal)
    Time.grid(row=3,column=2, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Время работы(τ с)')
    lb_t0.grid(row=2,column=2)
    r = Entry(Thermal)
    r.grid(row=5,column=2, padx=5, pady=5)
    lb_r = Label(Thermal, text='Сопротивление цепи (R Ом)')
    lb_r.grid(row=4,column=2)
    CapcTemp = Entry(Thermal)
    CapcTemp.grid(row=7,column=2, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Теплоемкость(с Дж/кг*°С)')
    lb_t0.grid(row=6,column=2)
    Btn_take = Button(Thermal, text="Расчёт", command = T25)
    Btn_take.grid(row=9,column=0)
    Power = Entry(Thermal)
    Power.grid(row=9,column=2, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Мощность(P Вт)')
    lb_t0.grid(row=8,column=2)
    '''Mt = Entry(Thermal)
    Mt.grid(row=10,column=0, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Температура плавления/парообразования/(t)')
    lb_t0.grid(row=9,column=0)'''
    Mass = Entry(Thermal)
    Mass.grid(row=9,column=1, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Масса(m кг)')
    lb_t0.grid(row=8,column=1)
    '''CapcTemp1 = Entry(Thermal)
    CapcTemp1.grid(row=12,column=1, padx=5, pady=5)
    lb_t0 = Label(Thermal, text='Теплоемкость после плавления/парообразования/(c)')
    lb_t0.grid(row=11,column=1)'''
    
    
def T25():
    global Temperature
    global Temperature0
    global Temperature1
    global Force
    global Power
    global Speed
    global Eff
    global Amperage
    global Voltage
    global Mass
    global CapcTemp
    global Mt
    global Time
    global CapcTemp1
    global r
    Mt = 999999
    PowerF = str(0)
    TemperatureF = str(0)
    QF = str(0)
    TimeF = str(0)
    rawPower = False
    UF = str(0)
    IF = str(0)
    RF = str(0)
    MassF = str(0)
    CapcTempF = str(0)
    if CapcTemp.get() != '':
        CapcTempF = float(CapcTemp.get())
    if Mass.get() != '':
        MassF = float(Mass.get())
    if Voltage.get() != '':
        UF = float(Voltage.get())
    if Amperage.get() != '':
        IF = float(Amperage.get())

    PowerF = Power.get()
    if Voltage.get() != '' and Amperage.get() != '':
        messagebox.showinfo("Сопротивление", str(float(Voltage.get())/float(Amperage.get()))+'Ом; R = U/I')
        RF = float(Voltage.get())/float(Amperage.get())
    if Voltage.get() != '' and r.get() != '':
        IF = float(Voltage.get())/float(r.get())
    if Amperage.get() != '' and r.get() != '':
        UF = float(Amperage.get())*float(r.get())
    if PowerF != '':
        PowerF = float(Power.get())
    if type(RF) != str and type(IF) != str:
        PowerF = RF**2*IF
    if type(RF) != str and type(UF) != str:
        PowerF = UF**2/RF
    if type(IF) != str and type(UF) != str:
        PowerF = float(UF*IF)
        messagebox.showinfo("Мощность тока", str(PowerF)+'Вт; P = U*I')
        rawPower = True
        if Power.get() != '':
            if float(Power.get()) > PowerF:
                messagebox.showinfo("КПД", str(100*PowerF/float(Power.get()))+'%; η = U*I/P')
            else:
                messagebox.showinfo("КПД", str(100*float(Power.get())/PowerF)+'%; η = P/U*I')
            
        
        
        
    if Temperature.get() == '' and Temperature1.get() != '' and Temperature0.get() != '':
        print(1)
        TemperatureF = float(Temperature1.get())-float(Temperature0.get())
        print(TemperatureF)
        messagebox.showinfo("Изменение температуры", str(-float(Temperature0.get())+float(Temperature1.get()))+'°С; dt = t1 - t0')
        k = (float(Temperature1.get())-float(Mt.get()))*(float(Temperature0.get())-float(Mt.get()))
    elif Temperature.get() != '':
        TemperatureF = float(Temperature.get())
        k = 1
    k = 1
    if type(TemperatureF) == float and  k > 0 and Mass.get() != '' and CapcTemp.get() != '':
        QF = TemperatureF*float(CapcTemp.get())*float(Mass.get())
        print(QF)
        messagebox.showinfo("Затраченное тепло", str(QF)+'Дж; Q = c*m*dt')
        if Time.get() != '':
            PowerF = QF/float(Time.get())
            messagebox.showinfo("Мощность", str(PowerF)+'Вт; Q/τ')
    if PowerF != '' and Time.get() != '':
        
        if Eff.get() != '':
            PowerF = PowerF/100*float(Eff.get())
        QF = PowerF*float(Time.get())
        messagebox.showinfo("Затраченное тепло", str(QF)+'Дж; Q = Q1/η')
        if Temperature0.get() == '' and Temperature1.get() != '':
            if k > 0 and CapcTemp.get() != '' and Mass.get() != '':
                TemperatureF = QF/(float(CapcTemp.get())*float(Mass.get()))
                messagebox.showinfo("Изменение температуры (°С)", TemperatureF)
                messagebox.showinfo("Начальная температура (°С)", float(Temperature1.get())-TemperatureF)
    elif PowerF != '' and type(QF) != str:
        
        if Eff.get() != '':
            PowerF = PowerF/100*float(Eff.get())
        TimeF = QF/PowerF
        messagebox.showinfo("Время работы", str(TimeF)+'c; τ = Q/P ')
        if Temperature0.get() == '' and Temperature1.get() != '':
            if k > 0 and CapcTemp.get() != '' and Mass.get() != '':
                TemperatureF = QF/(float(CapcTemp.get())*float(Mass.get()))
                messagebox.showinfo("Изменение температуры (°С)", TemperatureF)
                messagebox.showinfo("Начальная температура (°С)", float(Temperature1.get())-TemperatureF) 
                
            
    
    
    
    
    if Force.get() != '' and Speed.get() != '' and PowerF == '':
        PowerF = float(Force.get())*float(Speed.get())
        ForceF = float(Force.get()) 
        SpeedF = float(Speed.get())
        print(PowerF)
        messagebox.showinfo("Мощность", str(PowerF)+'Вт; P = F*V')
    elif Force.get() == '' and Speed.get() != '' and PowerF != '':
        ForceF = PowerF/float(Speed.get())
        SpeedF = float(Speed.get())
        print(ForceF)
        messagebox.showinfo("Сила Тяги", str(ForceF) + 'H; F = P/V')
    elif Force.get() != '' and Speed.get() == '' and PowerF != '':
        SpeedF = PowerF/float(Force.get())
        ForceF = float(Force.get())
        print(SpeedF)
        messagebox.showinfo("Скорость", str(SpeedF)+'м/с; V = P/F')
    elif Force.get() != '' and Speed.get() != '' and PowerF != '':
        PowerF1 = float(Force.get())*float(Speed.get())
        messagebox.showinfo("Мощность", str(PowerF1)+'Вт; P = F*V')
        if Power.get() != '':
            Power0 = float(Power.get())
        else:
            Power0 = PowerF
        if float(Power0) > PowerF1:
            messagebox.showinfo("КПД", str(100*PowerF1/Power0)+'%; η = F*V/P')
        else:
            messagebox.showinfo("КПД", str(100*Power0/PowerF1)+'%; η = P/F*V')
    if type(PowerF) == float:
        if Eff.get() != '' and rawPower == False:
            PowerF = PowerF*100/float(Eff.get())
            print(PowerF)
            messagebox.showinfo("Мощность до учета КПД", str(PowerF)+'Вт; P = P1/η')
        if Voltage.get() == '' and Amperage.get() != '':
            UF = PowerF / float(Amperage.get())
            print(UF)
            messagebox.showinfo("Напряжение", str(UF)+'В; U = F/I')
        elif Voltage.get() != '' and Amperage.get() == '':
            IF = PowerF / float(Voltage.get())
            print(IF)
            messagebox.showinfo("Сила тока", str(IF)+'А; I = F/U')
    if Voltage.get() != '' and Amperage.get() != '':
        PowerF = float(Voltage.get())*float(Amperage.get())
    if Mass.get() == '' and CapcTemp.get() != '' and type(QF) != str:
        MassF = QF/(float(CapcTemp.get())*TemperatureF)
        messagebox.showinfo("Масса", str(MassF)+'кг; m = Q/(c*dt)')
    if Time.get() != '':
        TimeF = float(Time.get())
    if type(TimeF) != str and type(PowerF) != str and type(TemperatureF) != str and CapcTemp.get() != '' and Mass.get() == '':
        messagebox.showinfo("Масса воды", str(TimeF*PowerF/(float(CapcTemp.get())*TemperatureF))+'кг; P*τ/c*dt')
        MassF =TimeF*PowerF/(float(CapcTemp.get()*TemperatureF))
    if type(UF) != str and type(IF) != str and r.get() == '' and type(RF) == str:
        messagebox.showinfo("Сопротивление", str(UF/IF)+'Ом; R = U/I')
    if type(PowerF) != str and type(TimeF) != str and type(MassF) != str and type(CapcTempF) != str:
        if Eff.get() == '' and type(TemperatureF) != str:
            messagebox.showinfo("КПД", str((TemperatureF*MassF*CapcTempF*100)/(PowerF*TimeF))+'%; η = 100%*c*m*dt/P*τ')
        elif Temperature.get() == '':
            TemperatureF = PowerF*TimeF/(MassF*CapcTempF)
            messagebox.showinfo("Разница температур", str(TemperatureF)+'°С; dt = P*τ/c*m')
    
    if Temperature1.get() == '' and Temperature0.get() != '' and type(TemperatureF) != str:
        messagebox.showinfo("Конечная темпераутра", str(TemperatureF+float(Temperature0.get()))+'°С; t1 = t0 + dt')
    elif Temperature1.get() != '' and Temperature0.get() == '' and type(TemperatureF) != str:
        messagebox.showinfo("Начальная темпераутра", str(-TemperatureF+float(Temperature1.get()))+'°С; t0 = t1 - dt')
    if type(PowerF) != str and type(TimeF) != str and type(MassF) != str and type(TemperatureF) != str and CapcTemp.get() == '':
        messagebox.showinfo("Теплоемкость", str((PowerF*TimeF)/(TemperatureF*MassF))+'Дж/кг*°С; c =  P*τ/dt*m')
        CapcTempF = (PowerF*TimeF)/(TemperatureF*MassF)
    if type(PowerF) != str and type(CapcTempF) != str and type(MassF) != str and type(TemperatureF) != str and Time.get() == '':
        messagebox.showinfo("Время", str(CapcTempF*TemperatureF*MassF/PowerF)+'с; τ =  с*dt*m/P')       
        
        
    
    
def Calc():
    global clicks
    global Sum_Hc
    global Energy
    global S_Mass
    global Melt_Data
    global Temperatur
    Max_Temp = Energy*Sum_Hc/sum(Mass)
    Melt_Data.sort()
    if Max_Temp < Melt_Data[0]:
        messagebox.showinfo("Ответ", Max_Temp)
    
    
def Clear():
    global clicks
    global Sum_Hc
    global Energy
    global S_Mass
    S_Mass.clear()
    S_Hc.clear()
    Energy = 0
    clicks = 0


def Activation():
    print(Energy, S_Hc, clicks, S_Mass)
    print(Melt_Data)


def Btn_Kin():
    global WindowCount
    WindowCount+=1
    if WindowCount > 1:
        tab_control.forget(1)
        WindowCount = 1
    Kin = ttk.Frame(tab_control)
    tab_control.add(Kin, text='Формулы')
    '''angle = Entry(Kin)
    angle.grid(row=0,column=2, padx=5, pady=5)
    lb_angle = Label(Kin, text='Угол')
    lb_angle.grid(row=0,column=1)
    frict_coeff = Entry(Ki', command = Calc_Kin)
    Btn_A.grid(column=0, row=0)
    Mass = Entry(Kin)
    Mass.grid(row=2,column=2, padx=5, pady=5)
    lb_m = Label(Kin, text='Масса')
    lb_m.grid(row=2,column=1)
    G = Entry(Kin)
    G.grid(row=3,column=2, padx=5, pady=5)
    lb_G = Label(Kin, text='Уск. свободного падения')
    lb_G.grid(row=3,column=1)
    N = math.cos(float(angle.get())/180*math.pi)*Massn)
    frict_coeff.grid(row=1,column=2, padx=5, pady=5)
    lb_frict_coeff = Label(Kin, text='Коэф.Трения')
    lb_frict_coeff.grid(row=1,column=1)
    Btn_A = Button(Kin, text='Расчет'''
    Q = Label(Kin, text='Q = c*m*dt = λ*m = L*m')
    Q.grid(row=0,column=1)
    Q = Label(Kin, text='A = P*t = Q/η')
    Q.grid(row=1,column=1)
    Q = Label(Kin, text='P = U*I = F*V')
    Q.grid(row=2,column=1)
    Q = Label(Kin, text='U = I*R')
    Q.grid(row=3,column=1)
    
    print(N)


def Calc_Kin():
    print('bob')
    

#Btn_Thermal = Button(Menu, text='Знт', command = Btn_Thermal)
#Btn_Thermal.grid(column=0, row=1)

Btn_Thermal = Button(Menu, text='Формулы', command = Btn_Kin)
Btn_Thermal.grid(column=1, row=0)

Btn_Thermal = Button(Menu, text='Калькулятор', command = Btn_С)
Btn_Thermal.grid(column=0, row=0)

#Btn_A = Button(Menu, text='Акт', command = Activation)
#Btn_A.grid(column=1, row=1)


tab_control.pack(expand=1, fill='both')


root.mainloop()