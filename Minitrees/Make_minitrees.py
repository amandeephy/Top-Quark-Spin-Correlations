import ROOT
import uproot
import argparse
import numpy as np
from   array import array
import matplotlib.pyplot as plt


def deltaphi(e_phi, m_phi):
    d_phi = e_phi - m_phi
    if (d_phi > np.pi):
        d_phi -= 2*np.pi
    if (d_phi < -np.pi):
        d_phi += 2*np.pi
    return d_phi


def dR(e_phi, e_eta, m_phi, m_eta):
    d_eta = abs(e_eta - m_eta)
    d_phi = deltaphi(e_phi, m_phi)
    return np.sqrt(d_phi**2 + d_eta**2)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input  Delphes Ntuple location')
    parser.add_argument('-o', '--output', help='Output Delphes Minitree location')

    args = parser.parse_args()
    inputFile = args.input
    outputFile = args.output

    fileptr = uproot.open(inputFile)['Delphes_Ntuples']

    # Jet MET
    jet_pt   = fileptr['jet_pt'].array()
    jet_eta  = fileptr['jet_eta'].array()
    jet_phi  = fileptr['jet_phi'].array()
    jet_mass = fileptr['jet_mass'].array()
    jet_btag = fileptr['jet_btag'].array()

    met_pt    = fileptr['met_pt'].array()
    met_phi   = fileptr['met_phi'].array()
    weight    = fileptr['weight'].array()
    scalar_ht = fileptr['scalar_ht'].array()

    # Electrons
    elec_pt     = fileptr['elec_pt'].array()
    elec_eta    = fileptr['elec_eta'].array()
    elec_phi    = fileptr['elec_phi'].array()
    elec_mass   = fileptr['elec_mass'].array()
    elec_charge = fileptr['elec_charge'].array()
    elec_reliso = fileptr['elec_reliso'].array()

    # Muons
    muon_pt     = fileptr['muon_pt'].array()
    muon_eta    = fileptr['muon_eta'].array()
    muon_phi    = fileptr['muon_phi'].array()
    muon_mass   = fileptr['muon_mass'].array()
    muon_charge = fileptr['muon_charge'].array()
    muon_reliso = fileptr['muon_reliso'].array()

    # Gen level jets
    genjet_pt   = fileptr['genjet_pt'].array()
    genjet_eta  = fileptr['genjet_eta'].array()
    genjet_phi  = fileptr['genjet_phi'].array()
    genjet_mass = fileptr['genjet_mass'].array()
    #genjet_btag = fileptr['genjet_btag'].array()

    # Gen level data
    genpart_pt     = fileptr['genpart_pt'].array()
    genpart_eta    = fileptr['genpart_eta'].array()
    genpart_phi    = fileptr['genpart_phi'].array()
    genpart_mass   = fileptr['genpart_mass'].array()
    genpart_pid    = fileptr['genpart_pid'].array()
    genpart_status = fileptr['genpart_status'].array()
    genpart_charge = fileptr['genpart_charge'].array()

    # Empty arrays
    # Similar to histograms in ROOT
 
    # By leading and subleading Pt
    l_pt   = []
    l_eta  = []
    l_phi  = []
    l_mass = []

    sl_pt   = []
    sl_eta  = []
    sl_phi  = []
    sl_mass = []

    # By Lepton flavor
    e_pt     = []
    e_eta    = []
    e_phi    = []
    e_charge = []

    mu_pt     = []
    mu_eta    = []
    mu_phi    = []
    mu_charge = []

    ljet_pt   = []
    ljet_eta  = []
    ljet_phi  = []
    ljet_mass = []

    sljet_pt   = []
    sljet_eta  = []
    sljet_phi  = []
    sljet_mass = []

    ST = []
    HT = []
    HT_check = []
    ET_miss  = []
    MET_phi  = []

    # Let's create a mask
    selection = np.zeros(len(jet_pt))

    # Loop over the events
    for i in range(len(jet_pt)):

        if (i % 1000 == 0):
            print('Processing event ' + str(i) + ' of ' + str(len(jet_pt)))

        # Temporary variables

        e_idx  = []
        mu_idx = []

        ef_idx  = []
        muf_idx = []

        jet_idx = []

        btag_cnt = 0

        e_4vec  = ROOT.TLorentzVector()
        mu_4vec = ROOT.TLorentzVector()

        ########## Electrons ##########

        # Ensure pt > 20 GeV and eta < 2.4 and isolation
        for j in range(len(elec_pt[i])):
            if (elec_pt[i][j] < 20):
                continue

            if ( abs(elec_eta[i][j]) > 2.4 or (1.4442 < abs(elec_eta[i][j]) < 1.5660)):
                continue

            #if (elec_reliso[i][j] > 0.0588):
                #continue

            e_idx.append(j)

        ###########  Muons ############

        # Ensure pt > 20 GeV and eta < 2.4 and isolation
        for j in range(len(muon_pt[i])):
            if (muon_pt[i][j] < 20):
                continue

            if (abs(muon_eta[i][j]) > 2.4):       # Should be 4?
                continue

            if (muon_reliso[i][j] > 0.15):
                continue

            mu_idx.append(j)

        # Ensure atleast one muon and one electron
        if (len(e_idx) == 0 or len(mu_idx) == 0):
            continue

        # Check for opp sign charge pairings
        for j in range(len(e_idx)):
            for k in range(len(mu_idx)):
                # e_idx and mu_idx have the list of valid electron and muon indices
                tmp_e_idx  = e_idx[j]
                tmp_mu_idx = mu_idx[k]

                if (elec_charge[i][tmp_e_idx] * muon_charge[i][tmp_mu_idx] == -1):
                    ef_idx.append(tmp_e_idx)
                    muf_idx.append(tmp_mu_idx)

        # Ensure such a pairing exists
        if (len(ef_idx) == 0 or len(muf_idx) == 0):
            continue

        # Assign leading indices to e and mu
        e_index = ef_idx[0]
        mu_index = muf_idx[0]

        # Defining the 4 vectors
        e_4vec.SetPtEtaPhiM(elec_pt[i][e_index]  , elec_eta[i][e_index] , elec_phi[i][e_index] , 0.000511)
        mu_4vec.SetPtEtaPhiM(muon_pt[i][mu_index], muon_eta[i][mu_index], muon_phi[i][mu_index], 0.105)

        # Mll cut (Step 3 according to the FW)
        if ((e_4vec + mu_4vec).M() < 20):
            continue

        ###########  Jets ###############

        for j in range(len(jet_pt[i])):

            # Ensure pt > 30 GeV and eta < 2.4 isolation

            if ((dR(elec_phi[i][e_index],  elec_eta[i][e_index], jet_phi[i][j], jet_eta[i][j]) < 0.4)
            or (dR(muon_phi[i][mu_index], muon_eta[i][mu_index], jet_phi[i][j], jet_eta[i][j]) < 0.4)):
                continue

            if ((jet_pt[i][j] < 30)):
                continue

            if ((abs(jet_eta[i][j]) > 2.4)):
                continue

            if (jet_btag[i][j] != 0):
                btag_cnt += 1

            jet_idx.append(j)

        # 2 Jets (Step 5 according to the FW)
        if(len(jet_idx) < 2):
            continue

        # Atleast one b-tag
        if (btag_cnt == 0):
            continue

        ljet_idx = jet_idx[0]
        sljet_idx = jet_idx[1]

        ####### Fill the arrays ##########

        # Leading and sub-leading lepton pts
        if (elec_pt[i][e_index] > muon_pt[i][mu_index] and elec_pt[i][e_index] > 25):
            l_pt.append(elec_pt[i][e_index])
            l_phi.append(elec_phi[i][e_index])
            l_eta.append(elec_eta[i][e_index])
            l_mass.append(elec_mass[i][e_index])

            sl_pt.append(muon_pt[i][mu_index])
            sl_phi.append(muon_phi[i][mu_index])
            sl_eta.append(muon_eta[i][mu_index])
            sl_mass.append(muon_mass[i][mu_index])

        elif (muon_pt[i][mu_index] > elec_pt[i][e_index] and muon_pt[i][mu_index] > 25):
            sl_pt.append(elec_pt[i][e_index])
            sl_phi.append(elec_phi[i][e_index])
            sl_eta.append(elec_eta[i][e_index])
            sl_mass.append(elec_mass[i][e_index])

            l_pt.append(muon_pt[i][mu_index])
            l_phi.append(muon_phi[i][mu_index])
            l_eta.append(muon_eta[i][mu_index])
            l_mass.append(muon_mass[i][mu_index])

        else:
            continue

        # By flavor
        e_pt.append(elec_pt[i][e_index])
        e_eta.append(elec_eta[i][e_index])
        e_phi.append(elec_phi[i][e_index])
        e_charge.append(elec_charge[i][e_index])

        mu_pt.append(muon_pt[i][mu_index])
        mu_eta.append(muon_eta[i][mu_index])
        mu_phi.append(muon_phi[i][mu_index])
        mu_charge.append(muon_charge[i][mu_index])

        # Leading and Subleading Pt
        ljet_pt.append(jet_pt[i][ljet_idx])
        ljet_phi.append(jet_phi[i][ljet_idx])
        ljet_eta.append(jet_eta[i][ljet_idx])
        ljet_mass.append(jet_mass[i][ljet_idx])

        sljet_pt.append(jet_pt[i][sljet_idx])
        sljet_phi.append(jet_phi[i][sljet_idx])
        sljet_eta.append(jet_eta[i][sljet_idx])
        sljet_mass.append(jet_mass[i][sljet_idx])

        ET_miss.append(met_pt[i][0])
        MET_phi.append(met_phi[i][0])

        # HT and ST
        # temporary variables
        Ht = 0
        St = 0

        for j in range(len(jet_pt[i])):
            Ht += jet_pt[i][j]
            St += jet_pt[i][j]

        for j in range(len(elec_pt[i])):
            St += elec_pt[i][j]

        for j in range(len(muon_pt[i])):
            St += muon_pt[i][j]

        St += met_pt[i][0]

        # Append to lists
        HT.append(Ht)
        ST.append(St)
        HT_check.append(scalar_ht[i][0])

        # Create a mask for selection
        selection[i] = 1

    # Other variables that we need
    llbar_deta = abs(np.array(l_eta) - np.array(sl_eta))
    llbar_dphi = abs(abs(abs(np.array(l_phi) - np.array(sl_phi)) - np.pi) - np.pi)

    bbbar_deta = abs(np.array(ljet_eta) - np.array(sljet_eta))
    bbbar_dphi = abs(abs(abs(np.array(ljet_phi) - np.array(sljet_phi)) - np.pi) - np.pi)

    # Store all info for additional jets and gen particles

    weight_sel = weight[selection == 1]

    jet_pt_sel   = jet_pt[selection == 1]
    jet_eta_sel  = jet_eta[selection == 1]
    jet_phi_sel  = jet_phi[selection == 1]
    jet_mass_sel = jet_mass[selection == 1]
    jet_btag_sel = jet_btag[selection == 1]

    gen_pt_sel     = genpart_pt[selection == 1]
    gen_pid_sel    = genpart_pid[selection == 1]
    gen_eta_sel    = genpart_eta[selection == 1]
    gen_phi_sel    = genpart_phi[selection == 1]
    gen_mass_sel   = genpart_mass[selection == 1]
    gen_status_sel = genpart_status[selection == 1]
    gen_charge_sel = genpart_charge[selection == 1]

    genjet_pt_sel   = genjet_pt[selection == 1]
    genjet_eta_sel  = genjet_eta[selection == 1]
    genjet_phi_sel  = genjet_phi[selection == 1]
    genjet_mass_sel = genjet_mass[selection == 1]
    #genjet_btag_sel = genjet_btag[selection == 1]

    # Make new root file with new tree
    opfile = ROOT.TFile(outputFile, 'recreate')
    tree   = ROOT.TTree("Step8", "Step8")
    hist   = ROOT.TH1F('Nevents', 'Nevents', 1, 0.5, 1.5)

    # Request arrays for branches

    maxn    = 9999
    HT_arr  = array('f', [0.])
    ST_arr  = array('f', [0.])
    MET_arr = array('f', [0.])
    MET_phi_arr  = array('f', [0.])
    HT_check_arr = array('f', [0.])

    l_pt_arr   = array('f', [0.])
    l_eta_arr  = array('f', [0.])
    l_phi_arr  = array('f', [0.])
    l_mass_arr = array('f', [0.])

    sl_pt_arr   = array('f', [0.])
    sl_eta_arr  = array('f', [0.])
    sl_phi_arr  = array('f', [0.])
    sl_mass_arr = array('f', [0.])

    # By flavor
    e_pt_arr      = array('f', [0.])
    e_eta_arr     = array('f', [0.])
    e_phi_arr     = array('f', [0.])
    e_charge_arr  = array('f', [0.])

    mu_pt_arr     = array('f', [0.])
    mu_eta_arr    = array('f', [0.])
    mu_phi_arr    = array('f', [0.])
    mu_charge_arr = array('f', [0.])

    ljet_pt_arr   = array('f', [0.])
    ljet_eta_arr  = array('f', [0.])
    ljet_phi_arr  = array('f', [0.])
    ljet_mass_arr = array('f', [0.])

    sljet_pt_arr   = array('f', [0.])
    sljet_eta_arr  = array('f', [0.])
    sljet_phi_arr  = array('f', [0.])
    sljet_mass_arr = array('f', [0.])

    llbar_deta_arr = array('f', [0.])
    llbar_dphi_arr = array('f', [0.])
    bbbar_deta_arr = array('f', [0.])
    bbbar_dphi_arr = array('f', [0.])

    weight_size_arr = array('i', [0])
    weight_arr      = array('f', maxn*[0.])

    jet_size_arr = array('i', [0])

    jet_pt_arr   = array('f', maxn*[0.])
    jet_eta_arr  = array('f', maxn*[0.])
    jet_phi_arr  = array('f', maxn*[0.])
    jet_mass_arr = array('f', maxn*[0.])
    jet_btag_arr = array('i', maxn*[0])

    genpart_size_arr = array('i', [0])

    genpart_pt_arr   = array('f', maxn*[0.])
    genpart_eta_arr  = array('f', maxn*[0.])
    genpart_phi_arr  = array('f', maxn*[0.])
    genpart_pid_arr  = array('i', maxn*[0])
    genpart_mass_arr = array('f', maxn*[0.])
    genpart_status_arr = array('i', maxn*[0])
    genpart_charge_arr = array('f', maxn*[0])

    genjet_size_arr = array('i', [0])

    genjet_pt_arr   = array('f', maxn*[0.])
    genjet_eta_arr  = array('f', maxn*[0.])
    genjet_phi_arr  = array('f', maxn*[0.])
    genjet_mass_arr = array('f', maxn*[0.])
    #genjet_btag_arr = array('i', maxn*[0])

    # Create the branches and assign the fill-variables to them as floats (F)

    tree.Branch("HT", HT_arr, 'HT/F')
    tree.Branch("ST", ST_arr, 'ST/F')
    tree.Branch("MET", MET_arr, 'MET/F')
    tree.Branch("HT_check", HT_check_arr, 'HT_check/F')
    tree.Branch("MET_phi" , MET_phi_arr , 'MET_phi/F')

    # Leading and sub-leading leptons
    tree.Branch("l_pt", l_pt_arr, 'l_pt/F')
    tree.Branch("l_eta", l_eta_arr, 'l_eta/F')
    tree.Branch("l_phi", l_phi_arr, 'l_phi/F')
    tree.Branch("l_mass", l_mass_arr, 'l_mass/F')

    tree.Branch("sl_pt", sl_pt_arr, 'sl_pt/F')
    tree.Branch("sl_eta", sl_eta_arr, 'sl_eta/F')
    tree.Branch("sl_phi", sl_phi_arr, 'sl_phi/F')
    tree.Branch("sl_mass", sl_mass_arr, 'sl_mass/F')

    # By flavor
    tree.Branch("e_pt"    , e_pt_arr    , 'e_pt/F')
    tree.Branch("e_eta"   , e_eta_arr   , 'e_eta/F')
    tree.Branch("e_phi"   , e_phi_arr   , 'e_phi/F')
    tree.Branch("e_charge", e_charge_arr, 'e_charge/F')

    tree.Branch("mu_pt"    , mu_pt_arr    , 'mu_pt/F')
    tree.Branch("mu_eta"   , mu_eta_arr   , 'mu_eta/F')
    tree.Branch("mu_phi"   , mu_phi_arr   , 'mu_phi/F')
    tree.Branch("mu_charge", mu_charge_arr, 'mu_charge/F')
 
    # Leading and sub-leading jets
    tree.Branch("ljet_pt", ljet_pt_arr, 'ljet_pt/F')
    tree.Branch("ljet_eta", ljet_eta_arr, 'ljet_eta/F')
    tree.Branch("ljet_phi", ljet_phi_arr, 'ljet_phi/F')
    tree.Branch("ljet_mass", ljet_mass_arr, 'ljet_mass/F')

    tree.Branch("sljet_pt", sljet_pt_arr, 'sljet_pt/F')
    tree.Branch("sljet_eta", sljet_eta_arr, 'sljet_eta/F')
    tree.Branch("sljet_phi", sljet_phi_arr, 'sljet_phi/F')
    tree.Branch("sljet_mass", sljet_mass_arr, 'sljet_mass/F')

    # Lepton angular stuff
    tree.Branch("llbar_deta", llbar_deta_arr, 'llbar_deta/F')
    tree.Branch("llbar_dphi", llbar_dphi_arr, 'llbar_dphi/F')
    tree.Branch("bbbar_deta", bbbar_deta_arr, 'bbbar_deta/F')
    tree.Branch("bbbar_dphi", bbbar_dphi_arr, 'bbbar_dphi/F')

    # Weights
    tree.Branch("weight_size", weight_size_arr, "weight_size/I")
    tree.Branch("weight", weight_arr, "weight[weight_size]/F")

    # Reco PUPPI jets
    tree.Branch("jet_size", jet_size_arr, "jet_size/I")
    tree.Branch("jet_btag", jet_btag_arr, "jet_btag[jet_size]/I")

    tree.Branch("jet_pt", jet_pt_arr, "jet_pt[jet_size]/F")
    tree.Branch("jet_eta", jet_eta_arr, "jet_eta[jet_size]/F")
    tree.Branch("jet_phi", jet_phi_arr, "jet_phi[jet_size]/F")
    tree.Branch("jet_mass", jet_mass_arr, "jet_mass[jet_size]/F")

    # Gen jets
    tree.Branch("genjet_size", genjet_size_arr, "genjet_size/I")
    #tree.Branch("genjet_btag", genjet_btag_arr, "genjet_btag[genjet_size]/I")

    tree.Branch("genjet_pt"  , genjet_pt_arr  , "genjet_pt[genjet_size]/F")
    tree.Branch("genjet_eta" , genjet_eta_arr , "genjet_eta[genjet_size]/F")
    tree.Branch("genjet_phi" , genjet_phi_arr , "genjet_phi[genjet_size]/F")
    tree.Branch("genjet_mass", genjet_mass_arr, "genjet_mass[genjet_size]/F")

    # Gen particles
    tree.Branch("genpart_size", genpart_size_arr, "genpart_size/I")
    tree.Branch("genpart_pid", genpart_pid_arr, "genpart_pid[genpart_size]/I")
    tree.Branch("genpart_status", genpart_status_arr,"genpart_status[genpart_size]/I")

    tree.Branch("genpart_pt", genpart_pt_arr, "genpart_pt[genpart_size]/F")
    tree.Branch("genpart_eta", genpart_eta_arr, "genpart_eta[genpart_size]/F")
    tree.Branch("genpart_phi", genpart_phi_arr, "genpart_phi[genpart_size]/F")
    tree.Branch("genpart_mass", genpart_mass_arr,"genpart_mass[genpart_size]/F")
    tree.Branch("genpart_charge", genpart_charge_arr,"genpart_charge[genpart_size]/F")

    for i in range(len(HT)):

        HT_arr[0]  = HT[i]
        ST_arr[0]  = ST[i]
        MET_arr[0] = ET_miss[i]

        HT_check_arr[0] = HT_check[i]
        MET_phi_arr[0]  = MET_phi[i]

        l_pt_arr[0]   = l_pt[i]
        l_eta_arr[0]  = l_eta[i]
        l_phi_arr[0]  = l_phi[i]
        l_mass_arr[0] = l_mass[i]

        sl_pt_arr[0]   = sl_pt[i]
        sl_eta_arr[0]  = sl_eta[i]
        sl_phi_arr[0]  = sl_phi[i]
        sl_mass_arr[0] = sl_mass[i]

        e_pt_arr[0]     = e_pt[i]
        e_eta_arr[0]    = e_eta[i]
        e_phi_arr[0]    = e_phi[i]
        e_charge_arr[0] = e_charge[i]

        mu_pt_arr[0]     = mu_pt[i]
        mu_eta_arr[0]    = mu_eta[i]
        mu_phi_arr[0]    = mu_phi[i]
        mu_charge_arr[0] = mu_charge[i]

        ljet_pt_arr[0]   = ljet_pt[i]
        ljet_eta_arr[0]  = ljet_eta[i]
        ljet_phi_arr[0]  = ljet_phi[i]
        ljet_mass_arr[0] = ljet_mass[i]

        sljet_pt_arr[0]   = sljet_pt[i]
        sljet_eta_arr[0]  = sljet_eta[i]
        sljet_phi_arr[0]  = sljet_phi[i]
        sljet_mass_arr[0] = sljet_mass[i]

        llbar_deta_arr[0] = llbar_deta[i]
        llbar_dphi_arr[0] = llbar_dphi[i]
        bbbar_deta_arr[0] = bbbar_deta[i]
        bbbar_dphi_arr[0] = bbbar_dphi[i]

        jet_size_arr[0]     = len(jet_pt_sel[i])
        weight_size_arr[0]  = len(weight_sel[i])
        genpart_size_arr[0] = len(gen_pt_sel[i])
        genjet_size_arr[0]  = len(genjet_pt_sel[i])

        for j in range(jet_size_arr[0]):
            jet_pt_arr[j]   = jet_pt_sel[i][j]
            jet_eta_arr[j]  = jet_eta_sel[i][j]
            jet_phi_arr[j]  = jet_phi_sel[i][j]
            jet_mass_arr[j] = jet_mass_sel[i][j]
            jet_btag_arr[j] = jet_btag_sel[i][j]

        for j in range(genjet_size_arr[0]):
            genjet_pt_arr[j]   = genjet_pt_sel[i][j]
            genjet_eta_arr[j]  = genjet_eta_sel[i][j]
            genjet_phi_arr[j]  = genjet_phi_sel[i][j]
            genjet_mass_arr[j] = genjet_mass_sel[i][j]
            #genjet_btag_arr[j] = genjet_btag_sel[i][j]

        for j in range(genpart_size_arr[0]):
            genpart_pt_arr[j]   = gen_pt_sel[i][j]
            genpart_pid_arr[j]  = gen_pid_sel[i][j]
            genpart_eta_arr[j]  = gen_eta_sel[i][j]
            genpart_phi_arr[j]  = gen_phi_sel[i][j]
            genpart_mass_arr[j] = gen_mass_sel[i][j]
            genpart_status_arr[j] = gen_status_sel[i][j]
            genpart_charge_arr[j] = gen_charge_sel[i][j]

        for k in range(weight_size_arr[0]):
            weight_arr[k] = weight_sel[i][k]

        tree.Fill()

    for i in range(len(selection)):
        hist.Fill(i)
    
    print('Number of events that pass selection :: ' + str(len(ljet_pt)))
    # Write the tree into the output file and close the file
    opfile.Write()
    opfile.Close()


main()
