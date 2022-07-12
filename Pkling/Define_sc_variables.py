import ROOT
import math
import uproot 
import argparse
import numpy  as np
import pandas as pd
from   datetime  import datetime

def ifOk(var_check) :
    if math.isfinite(var_check) :
        vOk = var_check
    else : 
        vOk = -999999.    
    return vOk


def main() :

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input SUSY root file')
    parser.add_argument('-o', '--output', help='Output SUSY pkl file')

    args       = parser.parse_args()
    inputFile  = args.input
    outputFile = args.output

    ptr      = uproot.open(inputFile)['Step8']  

    t_pt     = ptr['t_pt'].array()
    t_phi    = ptr['t_phi'].array()
    t_eta    = ptr['t_eta'].array()
    t_rap    = ptr['t_rapidity'].array()

    tbar_pt  = ptr['tbar_pt'].array()
    tbar_phi = ptr['tbar_phi'].array()
    tbar_eta = ptr['tbar_eta'].array()
    tbar_rap = ptr['tbar_rapidity'].array()

    tt_mass  = ptr['tt_mass'].array()

    l_pt     = ptr['l_pt'].array()
    l_phi    = ptr['l_phi'].array()
    l_eta    = ptr['l_eta'].array()
    l_mass   = ptr['l_mass'].array()

    lbar_pt   = ptr['lbar_pt'].array()
    lbar_phi  = ptr['lbar_phi'].array()
    lbar_eta  = ptr['lbar_eta'].array()
    lbar_mass = ptr['lbar_mass'].array()

    print('Finished loading in branches now to computing the spin corr')

    # The Bi's
    h_b1k = []
    h_b2k = []
    h_b1j = []
    h_b2j = []
    h_b1r = []
    h_b2r = []
    h_b1q = []
    h_b2q = []
    h_b1n = []
    h_b2n = []

    h_bP_kk = []
    h_bM_kk = []
    h_bP_jj = []
    h_bM_jj = []
    h_bP_rr = []
    h_bM_rr = []
    h_bP_qq = []
    h_bM_qq = []
    h_bP_nn = []
    h_bM_nn = []

    # The Ci's
    h_ckk = []
    h_crr = []
    h_cnn = []
    h_crk = []
    h_ckr = []
    h_cnr = []
    h_crn = []
    h_cnk = []
    h_ckn = []

    h_cP_rk = []
    h_cM_rk = []
    h_cP_nr = []
    h_cM_nr = []
    h_cP_nk = []
    h_cM_nk = []

    # ll variables
    h_ll_dphi = []
    h_ll_deta = []
    h_ll_dR   = []
    h_c_hel   = []

    # Adapted from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/Configuration/analysis/diLeptonic/src/VariablesPhiTT.cc

    for i in range(len(t_pt)) :

        if (i% 1000 == 0):
            print('Processing event :: ' + str(i))
            now          = datetime.now()   # Time keeping
            current_time = now.strftime("%H:%M:%S")
            print("Current Time =", current_time)
        
        top    = ROOT.TLorentzVector()
        atop   = ROOT.TLorentzVector() 
        lep    = ROOT.TLorentzVector()
        alep   = ROOT.TLorentzVector()
        
        top.SetPtEtaPhiM(t_pt[i]    , t_eta[i]   , t_phi[i]   , 172.5)
        atop.SetPtEtaPhiM(tbar_pt[i], tbar_eta[i], tbar_phi[i], 172.5)  
        
        lep.SetPtEtaPhiM(l_pt[i]    , l_eta[i]   , l_phi[i]   , l_mass[i])
        alep.SetPtEtaPhiM(lbar_pt[i], lbar_eta[i], lbar_phi[i], lbar_mass[i])

        # Spin corr variables without boosting
        ll_dEta = ifOk( lep.Eta() - alep.Eta() )
        ll_dPhi = ifOk( lep.DeltaPhi( alep ) )
        ll_dR   = ifOk( lep.DeltaR( alep ) )
        cHel    = ifOk( alep.Vect().Unit().Dot( lep.Vect().Unit() ) )

        # The various Bernreuther bases
        kBase = ROOT.TVector3()
        jBase = ROOT.TVector3()
        qBase = ROOT.TVector3()
        #rBase = ROOT.TVector3()
        #nBase = ROOT.TVector3()
        
        # Beam unit vector
        p3_pPro = ROOT.TVector3(0., 0., 1.)
        
        # The bases definition: Bernreuther spinMatrix 1508.05271
        p4_TT   = top + atop
        b4_TT   = ROOT.TVector3(-1. * p4_TT.BoostVector() )
        
        b4_pTop = top
        b4_pTop.Boost( b4_TT )
        
        b4_aTop = atop
        b4_aTop.Boost( b4_TT )
        
        # Maybe need to be careful with the signs here ?? Check how the pairings are implemented in the FW
        b4_aLep = alep
        b4_aLep.Boost( b4_TT )
        b4_aLep.Boost( -1. * b4_pTop.BoostVector() )
        
        b4_pLep = lep
        b4_pLep.Boost( b4_TT )
        b4_pLep.Boost( -1. * b4_aTop.BoostVector() )
        
        # Calculating the top-beam angle for pTop only
        c_pTP = b4_pTop.Vect().Unit().Dot(p3_pPro)
        s_pTP = np.sqrt(1. - (c_pTP * c_pTP))

        # The signs needed to account for Bose symmetry
        sY = 1. if ( c_pTP >= 0.) else -1.
        sD = 1. if ( abs(top.Rapidity()) >= abs(atop.Rapidity()) ) else -1. 

        # Define the base vectors a
        # j and q base are the k* and r* respectively
        # b is always -a
        
        kBase = b4_pTop.Vect().Unit()
        jBase = sD * kBase
        r_arr = (sY / s_pTP) * (p3_pPro - (c_pTP * kBase)) # Store in a temp np array since pyROOT typecasts automatically
        rBase = ROOT.TVector3(r_arr[0], r_arr[1], r_arr[2]).Unit()
        qBase = sD * rBase
        n_arr = (sY / s_pTP) *  p3_pPro.Cross(kBase)       # Store in a temp np array since pyROOT typecasts automatically
        nBase = ROOT.TVector3(n_arr[0], n_arr[1], n_arr[2]).Unit()

        # Find the relevant angles in these bases
        ck_aLep = b4_aLep.Vect().Unit().Dot( kBase )
        ck_pLep = b4_pLep.Vect().Unit().Dot( -1. * kBase )

        cj_aLep = b4_aLep.Vect().Unit().Dot( jBase )
        cj_pLep = b4_pLep.Vect().Unit().Dot( -1. * jBase )

        cr_aLep = b4_aLep.Vect().Unit().Dot( rBase )
        cr_pLep = b4_pLep.Vect().Unit().Dot( -1. * rBase )

        cq_aLep = b4_aLep.Vect().Unit().Dot( qBase )
        cq_pLep = b4_pLep.Vect().Unit().Dot( -1. * qBase )

        cn_aLep = b4_aLep.Vect().Unit().Dot( nBase )
        cn_pLep = b4_pLep.Vect().Unit().Dot( -1. * nBase )
        
        # Fill the raw angles into VarFloats
        b1k = ifOk( ck_aLep )
        b2k = ifOk( ck_pLep )

        b1j = ifOk( cj_aLep )
        b2j = ifOk( cj_pLep )

        b1r = ifOk( cr_aLep )
        b2r = ifOk( cr_pLep )

        b1q = ifOk( cq_aLep )
        b2q = ifOk( cq_pLep )

        b1n = ifOk( cn_aLep )
        b2n = ifOk( cn_pLep )

        # Now we can squeeze it all out based on table 5 page 16
        # The B1 ~ c_aLep, B2 ~ c_pLep sums
        
        bP_kk = ifOk( ck_aLep + ck_pLep )
        bM_kk = ifOk( ck_aLep - ck_pLep )

        bP_jj = ifOk( cj_aLep + cj_pLep )
        bM_jj = ifOk( cj_aLep - cj_pLep )

        bP_rr = ifOk( cr_aLep + cr_pLep )
        bM_rr = ifOk( cr_aLep - cr_pLep )

        bP_qq = ifOk( cq_aLep + cq_pLep )
        bM_qq = ifOk( cq_aLep - cq_pLep )

        bP_nn = ifOk( cn_aLep + cn_pLep )
        bM_nn = ifOk( cn_aLep - cn_pLep )
        
        # spinCorr coeff Cab = -9<cab>
        ckk = ifOk( ck_aLep * ck_pLep )
        crr = ifOk( cr_aLep * cr_pLep )
        cnn = ifOk( cn_aLep * cn_pLep )

        crk = cr_aLep * ck_pLep
        ckr = ck_aLep * cr_pLep

        cnr = cn_aLep * cr_pLep
        crn = cr_aLep * cn_pLep

        cnk = cn_aLep * ck_pLep
        ckn = ck_aLep * cn_pLep

        cP_rk = ifOk( crk + ckr )
        cM_rk = ifOk( crk - ckr )

        cP_nr = ifOk( cnr + crn )
        cM_nr = ifOk( cnr - crn )

        cP_nk = ifOk( cnk + ckn )
        cM_nk = ifOk( cnk - ckn )
        
        # Fill the empty lists
        h_b1k.append(b1k)
        h_b2k.append(b2k)
        h_b1j.append(b1j)
        h_b2j.append(b2j)
        h_b1r.append(b1r)
        h_b2r.append(b2r)
        h_b1q.append(b1q)
        h_b2q.append(b2q)
        h_b1n.append(b1n)
        h_b2n.append(b2n)

        h_bP_kk.append(bP_kk)
        h_bM_kk.append(bM_kk)
        h_bP_jj.append(bP_jj)
        h_bM_jj.append(bM_jj)
        h_bP_rr.append(bP_rr)
        h_bM_rr.append(bM_rr)
        h_bP_qq.append(bP_qq)
        h_bM_qq.append(bM_qq)
        h_bP_nn.append(bP_nn)
        h_bM_nn.append(bM_nn)
        
        h_ckk.append(ckk)
        h_crr.append(crr)
        h_cnn.append(cnn)
        h_crk.append(crk)
        h_ckr.append(ckr)
        h_cnr.append(cnr)
        h_crn.append(crn)
        h_cnk.append(cnk)
        h_ckn.append(ckn)
        
        h_cP_rk.append(cP_rk)
        h_cM_rk.append(cM_rk)
        h_cP_nr.append(cP_nr)
        h_cM_nr.append(cM_nr)
        h_cP_nk.append(cP_nk)
        h_cM_nk.append(cM_nk)
        
        h_ll_dphi.append(ll_dPhi)
        h_ll_deta.append(ll_dEta)
        h_ll_dR.append(ll_dR)
        h_c_hel.append(cHel)

    print('Storing to dataframe')

    arr        = [h_ckk  ,   h_ckr,   h_crk,   h_cnn,   h_cnk,   h_ckn,   h_crr,   h_crn,   h_cnr, 
                  h_cP_rk, h_cM_rk, h_cP_nr, h_cM_rk, h_cP_nk, h_cM_nk,
                  h_b1k  ,   h_b2k,   h_b1q,   h_b2q,   h_b1r,   h_b2r,   h_b1n,   h_b2n,   h_b1j,   h_b2j,     
                  h_bP_kk, h_bM_kk, h_bP_jj, h_bM_jj, h_bP_rr, h_bM_rr, h_bP_qq, h_bM_qq, h_bP_nn, h_bM_nn,  
                  h_ll_dphi,      h_ll_deta, h_ll_dR, h_c_hel]

    df         = pd.DataFrame(data = arr).T

    df.columns = [ 'c_kk',  'c_kr',  'c_rk',  'c_nn',  'c_nk',  'c_kn',  'c_rr',  'c_rn',  'c_nr', 
                    'cP_rk', 'cM_rk', 'cP_nr', 'cM_rk', 'cP_nk', 'cM_nk',
                    'b1k',   'b2k',   'b1q',   'b2q',   'b1r',   'b2r',   'b1n',   'b2n',   'b1j',   'b2j', 
                    'bP_kk', 'bM_kk', 'bP_jj', 'bM_jj', 'bP_rr', 'bM_rr', 'bP_qq', 'bM_qq', 'bP_nn', 'bM_nn',
                    'llbar_dphi', 'llbar_deta',    'dR', 'c_hel']


    df.to_pickle('Pkl_ips/SUSY_top_reco_op/' + str(outputFile))

    print('Done processing file')

main()
