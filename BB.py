import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(layout="wide")

st.write("""
## Bolig budsjett
""")


# Input

st.sidebar.header("Bolig og finnansiering")
# Bolig pris 
P_0 = st.sidebar.number_input('Boligpris',value=5500000,step=50000)
# Lånemengde
L_0 = st.sidebar.number_input('Lånebeløp',value=4000000,step=50000)
# Egenkapital
EK = st.sidebar.number_input('Egenkaptial',value=1500000,step=50000)
# Årslønn
lønn_0 = st.sidebar.number_input('Årslønn før skatt',value=550000,step=10000)
# Lengde på lån
lt = st.sidebar.selectbox('Låneperode',np.arange(30)+1, 24)

# Inndeling
st.sidebar.subheader("Detaljer")
    
# Rente
r = st.sidebar.slider('Rente i %', -5.0, 15.0, 3.0,step=0.1)
r = r/100
# Prisvekst bolig
p = st.sidebar.slider('Prisvekst bolig', -5.0, 15.0, 3.0,step=0.1)
p = p/100
# Markedsavkastning
m = st.sidebar.slider('Markedsavkastning', -5.0, 15.0, 5.0,step=0.1)
m = m/100
# KPI inflasjon
cpi = st.sidebar.slider('KPI inflasjon', -5.0, 15.0, 2.0,step=0.1)
cpi = cpi/100
# Årlig lønnsvekst
lg = st.sidebar.slider('Årlig lønnsvekst', -5.0, 15.0, 3.0,step=0.1)
lg = lg/100

# Inndeling
st.sidebar.subheader("Utgifter og inntekter")

# Andre Utgifter
au_kr = st.sidebar.slider('Felleskonstander +', 0, 10000, 5000,step=100)
# Månedlig utgifter
kon_0 = st.sidebar.slider('Månedsforbruk', 0, 25000, 13000,step=500)
# Leieinntekter
leitak_kr = st.sidebar.slider('Leieinntekter', 0, 25000, 6000,step=100)
# Studiegjeld 
sl = st.sidebar.slider('Studiegjeld nedbetaling per måned', 0, 10000, 3500,step=100)
# Tid til betjening av studiegjeld
t_sl = st.sidebar.slider('Studiegjeld utsettelse, måneder', 0, 36, 36)

# Inndeling
st.sidebar.subheader("Graf")
# Tidshorisont
t = st.sidebar.slider('Tidshorisont', 1, 30, 8)


# Programmer
def pris_S(P_0,p,lt):
    P_S = np.zeros((lt)*12+1)
    for i in np.arange((lt)*12+1):
        P_S[i] = P_0*(1+p/12)**(i)
    return P_S

def Salg_pris_S(P_S,salg_kost,cpi):
    # pris etter salgskostnader
    SK_S = np.zeros(len(P_S))
    for i in np.arange(len(P_S)):
        SK_S[i] = salg_kost*(1+cpi/12)**(i)
    Salg_pris_S = P_S-SK_S
    return Salg_pris_S

def lån_S_L(L_0,r,lt):
    # Størelse på Lån, serie
    L_S = np.zeros((lt)*12+1)
    A = ann(L_0,r,lt)
    for i in np.arange(lt*12+1):
        L_S[i] = lån_0(A,r,(lt)-i/12)
    return L_S

def ann(L_0,r,lt=25):
    # Annuitetsbetaling, mnd lånekostander
    A = L_0*(r/12*(1+r/12)**(lt*12))/((1+r/12)**(lt*12)-1)
    return A

def ann_S_L(L_0,r,lt=25):
    # annuitet serie med lån som innput
    A = ann(L_0,r,lt)
    A_S = ann_S_A(A,lt)
    A_S = np.insert(A_S,0,0)
    return A_S

def avd_S(L_S):
    # Avdrag på lån hver mnd
    avd_S = [L_S[i - 1] - x  for i, x in enumerate(L_S)][1:]
    avd_S.insert(0,0)
    return avd_S

def lønn_S(lønn_0,lg,lt):
    lønn_S = np.zeros(lt)
    for i in np.arange(lt):
        lønn_S[i] = lønn_0*(1+lg)**i
    lønn_S = (lønn_S/12).repeat(12) 
    lønn_S = np.insert(lønn_S,0,0)
    return lønn_S

def lønn_e_skatt_S(lønn_S):
    lønn_skatt_S = np.zeros(len(lønn_S))
    for i in np.arange(len(lønn_S)):
        lønn_skatt_S[i] = (lønn_S[i]*12*(0.1847+0.0000001396*lønn_S[i]*12))/12
    lønn_e_skatt_S = lønn_S-lønn_skatt_S
    return lønn_e_skatt_S

def kon_S(kon_0,cpi,lt):
    kon_S = np.zeros(lt*12)
    for i in np.arange(lt*12):
        kon_S[i] = kon_0*(1+cpi/12)**i
    kon_S = np.insert(kon_S,0,0)
    return kon_S

def ek_inv(inv_S,m):
    ek_inv_S = np.zeros(len(inv_S))
    for i in np.arange(len(inv_S)):
        ek_inv_S[i] = ek_inv_S[i-1]*(1+m/12)+inv_S[i]
    return ek_inv_S

def lån_0(A,r,lt=25):
    # Orginal låmnemengde
    L = A*((1+r/12)**(lt*12)-1)/(r/12*(1+r/12)**(lt*12))
    return L

def lån_S_A(A,r,lt):
    # Annutiet serie ????? 
    L_S = np.zeros((lt)*12+1)
    for i in np.arange((lt)*12+1):
        L_S[i] = lån(A,r,((lt)-i/12))
    return L_S

def ann_S_A(A, lt):
    # Annuitet serie med annuitet som innput 
    A_S = np.zeros(lt*12)
    for i in np.arange(lt*12):
        A_S[i] = A
    return A_S


def summary(P_0, L_0, EK, lønn_0, p, r, m, cpi, lt, kon_0, t, sl,lg,leitak_kr, au_kr,t_sl):
    verdi_sats="Sats"
    kap_skatt=0.27
    salg_kost=100000
    au=0.015
    leie_in=0
    # EK bolig, P_0, p, lt, L_0,r
    P_S = pris_S(P_0,p,lt)
    P_S_S = Salg_pris_S(P_S,salg_kost,cpi)
    L_S = lån_S_L(L_0,r,lt)
    EK_S = P_S_S-L_S # riktig

    # Utgifter
    A = ann(L_0,r,lt)
    A_S = ann_S_L(L_0,r,lt)

    AU_0 = au_kr
    AU_S= np.zeros((lt)*12+1)
    for i in np.arange((lt)*12+1):
        AU_S[i] = AU_0*(1+cpi/12)**(i)
    AU_S[0] = 0

    av_S = avd_S(L_S)
    r_k_S = A_S - av_S

    su_S = r_k_S+AU_S


    Studie_gjeld_S = np.concatenate(((np.repeat(0,t_sl+1)),(np.repeat(sl,20*12+1)),(np.repeat(0,(120-t_sl)+2))))
    Studie_gjeld_S = Studie_gjeld_S[:lt*12+1]


    # lønn, lønn_0,cpi,lt
    # lønn
    lønn_S_1 = lønn_S(lønn_0,lg,lt) # riktig
    # lønn etter skatt
    lønn_e_skatt_S_1 = lønn_e_skatt_S(lønn_S_1)
    # etter utgifter
    lønn_e_ut_1 = lønn_e_skatt_S_1 - AU_S - A_S - Studie_gjeld_S

    # leie, A, P_0, au, leie_in, cpi, lt
    leie_in_S_1 = np.zeros(lt)
    for i in np.arange(lt):
        leie_in_S_1[i] = (leitak_kr*(1+cpi)**i)
    leie_in_S_1 = leie_in_S_1.repeat(12)
    leie_in_S_1 = np.insert(leie_in_S_1,0,0)
        
        
    # Konsum, kon_0, cpi, lt
    kon_S_1 = kon_S(kon_0,cpi,lt) # riktig

    # investering, EK, m
    inv_S_1 = lønn_e_ut_1 + leie_in_S_1 - kon_S_1
    inv_S_1[0] = EK-EK_S[0]-salg_kost
    ek_inv_S = ek_inv(inv_S_1,m)
    ek_inv_es_S = ek_inv_S-((ek_inv_S-np.cumsum(inv_S_1))*kap_skatt)   # etter skatt

    sum_EK = EK_S+ek_inv_es_S
    
    # Format
    P_S = [ '%.0f' % elem for elem in P_S ]
    L_S = [ '%.0f' % elem for elem in L_S ]
    ek_inv_S = [ '%.0f' % elem for elem in ek_inv_S ]
    EK_S = [ '%.0f' % elem for elem in EK_S ]
    ek_inv_es_S = [ '%.0f' % elem for elem in ek_inv_es_S ]
    sum_EK = [ '%.0f' % elem for elem in sum_EK ]
    lønn_S_1 = [ '%.0f' % elem for elem in lønn_S_1 ]
    lønn_e_skatt_S_1 = [ '%.0f' % elem for elem in lønn_e_skatt_S_1 ]
    A_S = [ '%.0f' % elem for elem in A_S ]
    av_S = [ '%.0f' % elem for elem in av_S ]
    r_k_S = [ '%.0f' % elem for elem in r_k_S ]
    AU_S = [ '%.0f' % elem for elem in AU_S ]
    su_S = [ '%.0f' % elem for elem in su_S ]
    Studie_gjeld_S = [ '%.0f' % elem for elem in Studie_gjeld_S ]
    lønn_e_ut_1 = [ '%.0f' % elem for elem in lønn_e_ut_1 ]
    leie_in_S_1 = [ '%.0f' % elem for elem in leie_in_S_1 ]
    kon_S_1 = [ '%.0f' % elem for elem in kon_S_1 ]
    inv_S_1 = [ '%.0f' % elem for elem in inv_S_1 ]
    
    
    # oversikt combined
    col = ["Bolig verdi","Lån","Verdipapirer","EK Bolig","EK verdipapirer","EK total","Lønn","Lønn etter skatt","Betjening av lån","Avdrag","Renteutgivter","Andre utgifter", "Sum Utgifter","Studie gjeld","lønn etter utgift","Leieinntekter","Konsum","Sparing"]
    data =  np.transpose([P_S,L_S,ek_inv_S,EK_S,ek_inv_es_S,sum_EK,lønn_S_1,lønn_e_skatt_S_1,A_S,av_S,r_k_S,AU_S,su_S,Studie_gjeld_S,lønn_e_ut_1, leie_in_S_1,kon_S_1,inv_S_1])

    summar_all =  pd.DataFrame(data=data, columns=col)


    results = summar_all

        # begrenset område
    results = results[:(t*12+1)]

    return results
 

summary_dataframe = summary(P_0, L_0, EK, lønn_0, p, r, m, cpi, lt, kon_0, t, sl,lg,leitak_kr, au_kr, t_sl)





col_l = ["Sparing","Bolig verdi","Lån","Verdipapirer","EK Bolig","EK verdipapirer","EK total","Lønn","Lønn etter skatt","Betjening av lån","Andre utgifter","Studie gjeld","lønn etter utgift","Leieinntekter","Konsum"]

select_verdi = st.selectbox('Velg verdi',col_l)
st.write(select_verdi)


plot_data = (summary_dataframe[select_verdi].tolist())[1:]
plot_data = np.float_(plot_data)



st.line_chart(plot_data)

st.write(summary_dataframe)



