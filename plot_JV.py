import pandas as pd
import matplotlib.pyplot as plt
import simudo.example.fourlayer.sweep_extraction as se

def main():
    date = 'Mar31'
    filename = f"example_output/{date}/a/HMA_V.csv"
    df = pd.read_csv(filename, skiprows=[1]) #skip the row containing units

    df['V'] = df['sweep_parameter:V']
    df['J'] = df['avg:current_CB:p_contact'] + df['avg:current_VB:p_contact']

    fig, ax = plt.subplots()
    ax.plot(df['V'], -df['J'], label="$J(V)$, (start of project)", marker='.')
    ax.set_xlabel(r'applied potential ($\mathrm{V}$)')
    ax.set_ylabel(r'total current ($\mathrm{mA}/\mathrm{cm}^2$)')
    ax.legend()
    plt.ylim([0, -df['J'][0]+10])
    plt.grid()
    plt.savefig(f'{date}_JV.png')
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(df['V'], -df['J']*df['V'], label="$P(V)$", marker='.')
    #ax.set_yscale('log')
    ax.set_title(f'Without contact layer')
    ax.set_xlabel(r'applied potential ($\mathrm{V}$)')
    ax.set_ylabel(r'Output power ($\mathrm{W}/\mathrm{cm}^2$)')
    ax.legend()
    plt.ylim([0, max(-df['J']*df['V'])+5])
    plt.grid()
    plt.savefig(f'{date}_PV.png')
    plt.show()


    print('Jsc =', df['J'][0])
    voc_index = findroot(df['J'])

    Voc_bound1 = df['V'][voc_index]
    Voc_bound2 = df['V'][voc_index+1]
    print('Voc lower bound:', round(Voc_bound1,4), 'V')
    print('Voc upper bound:', round(Voc_bound2,4), 'V')
    print('Voc ~~', round((Voc_bound1+Voc_bound2)/2,4), 'V')


def findroot(lst):
    for i,j in enumerate(lst):
        if lst[i]*lst[i+1] < 0:
            return i
        if i+2 >= len(lst):
            raise Exception('###### ERROR: Voc NOT FOUND ######')
            return -1

if __name__ == '__main__':
    main()
