alt = float(input("What is the altitude in feet?\n"))

temp = 0

if alt <= 36089:
    temp = 518.69-3.5662e-3*alt
if 36089 < alt <= 65617:
    temp = 389.99
if alt > 65617:
    temp = 389.99+5.4864e-4*(alt-65617)

print("The Air Temperature is: ", temp)
