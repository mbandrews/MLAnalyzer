from __future__ import print_function
import numpy as np
np.random.seed(0)
import os, glob
import time
import argparse
import pickle
#from data_utils import *
import json
from collections import OrderedDict

#'''
import torch
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import *
import torch.nn as nn
#os.environ["CUDA_VISIBLE_DEVICES"]=str(0)
#'''

import ROOT
import pileupParser

# Register command line options
parser = argparse.ArgumentParser(description='Run STEALTH selection.')
#parser.add_argument('-s', '--sample', default=None, type=str, help='Sample name.')
parser.add_argument('-m', '--minbXsec', default=69200, type=float, help='min bias xsec [microbarn]. run2 default: 69200')
parser.add_argument('-p', '--inputLumiJSON', default='Collisions17_13TeV_PileUp_pileup_latest.txt', type=str, help='Input PU json file.')
parser.add_argument('-i', '--img_inputs', default='img_inputs.txt', type=str, help='List file of img inputs.')
#parser.add_argument('-i', '--inputs', default=['MAntuples/Run2017F_mantuple.root'], nargs='+', type=str, help='Input MA ntuple.')
parser.add_argument('-o', '--outdir', default='.', type=str, help='Output directory.')
args = parser.parse_args()

#############################
# Crop out EB shower from full EB image
def crop_EBshower(imgEB, ieta, iphi, window=32):

    # NOTE: image window here should correspond to the one used in RHAnalyzer
    off = window//2
    ieta = int(ieta[0])+1 # seed positioned at [15,15]
    iphi = int(iphi[0])+1 # seed positioned at [15,15]

    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,-diff:],
                                   imgEB[:,ieta-off:ieta+off,:iphi+off]), axis=-1)
    # Wrap-around on right side
    elif 360-iphi < off:
        diff = off - (360-iphi)
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,iphi-off:],
                                   imgEB[:,ieta-off:ieta+off,:diff]), axis=-1)
    # Nominal case
    else:
        img_crop = imgEB[:,ieta-off:ieta+off,iphi-off:iphi+off]

    return img_crop

# Crop out EB shower from full EB image with padding
def crop_EBshower_padded(imgEB, ieta, iphi, window=32):

    assert len(imgEB.shape) == 3, '!! len(imgEB.shape): %d != 3'%len(imgEB.shape)
    assert ieta[0] < imgEB.shape[1], '!! ieta:%d !< imgEB.shape[1]:%d'%(ieta[0], imgEB.shape[1])
    assert iphi[0] < imgEB.shape[2], '!! iphi:%d !< imgEB.shape[2]:%d'%(iphi[0], imgEB.shape[2])

    # NOTE: image window here should correspond to the one used in RHAnalyzer
    off = window//2
    ieta = int(ieta[0])+1 # seed positioned at [15,15]
    iphi = int(iphi[0])+1 # seed positioned at [15,15]

    # ------------------------------------------------
    # ieta (row) padding
    # ------------------------------------------------
    pad_lo, pad_hi = 0, 0
    # lower padding check
    if ieta >= off:
        ieta_lo = ieta-off
    else:
        pad_lo = abs(ieta-off)
        ieta_lo = 0
    # upper padding check
    if ieta+off <= imgEB.shape[1]:
        ieta_hi = ieta+off
    else:
        pad_hi = abs(ieta+off-imgEB.shape[1])
        ieta_hi = imgEB.shape[1]

    # ------------------------------------------------
    # iphi (col) wrap-around
    # ------------------------------------------------
    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop_ = np.concatenate((imgEB[:, ieta_lo:ieta_hi, -diff:],
                                    imgEB[:, ieta_lo:ieta_hi, :iphi+off]), axis=-1)
    # Wrap-around on right side
    elif 360-iphi < off:
        diff = off - (360-iphi)
        img_crop_ = np.concatenate((imgEB[:, ieta_lo:ieta_hi, iphi-off:],
                                    imgEB[:, ieta_lo:ieta_hi, :diff]), axis=-1)
    # Nominal case
    else:
        img_crop_ = imgEB[:, ieta_lo:ieta_hi, iphi-off:iphi+off]

    # Add ieta padding if needed
    img_crop = np.pad(img_crop_, ((0,0), (pad_lo, pad_hi), (0,0)), 'constant') # pads with 0
    assert img_crop.shape[1] == window, '!! img_crop.shape[1]:%d != window:%d'%(img_crop.shape[1], window)
    assert img_crop.shape[2] == window, '!! img_crop.shape[2]:%d != window:%d'%(img_crop.shape[2], window)

    return img_crop

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

#sample = args.sample
sample = args.img_inputs.split('_')[-1].split('.txt')[0]
print('>> Sample:', sample)
#sample = 'DoublePi0Pt10To100_m0To1600_pythia8_ReAOD_PU2017_MINIAODSIM_pu'

#############################
# PU info
minbXsec = args.minbXsec # microbarn, run2 default: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
inputLumiJSON = args.inputLumiJSON

def parseInputFile(inputfilename):
    '''
    output ({run:[ls:[inlumi, meanint]]})
    '''
    selectf = open(inputfilename,'r')
    inputfilecontent = selectf.read()
    p = pileupParser.pileupParser(inputfilecontent)

    runlsbyfile = p.runsandls()
    return runlsbyfile

def get_pu(runId, lumiId, inputPileupRange):

    pu = -1.
    if runId in inputPileupRange.keys():
        if lumiId in inputPileupRange[runId]:
            lumiInfo = inputPileupRange[runId][lumiId]
            pu = lumiInfo[2]*minbXsec

    return pu

inputPileupRange = parseInputFile(inputLumiJSON)

# for MC only
def get_pumc(tree):

    bx = np.array(tree.puBX)
    pu = np.array(tree.puTrue)
    ibx0 = np.argwhere(bx == 0)
    return pu[ibx0]

#############################
# Load IMG ntuples as main TTree
'''
eos_basedir = '/store/user/lpcml/mandrews/IMG'
img_inputs = run_eosfind(eos_basedir, sample)
'''
print('Opening img input list:',args.img_inputs)
img_inputs = []
with open(args.img_inputs, 'r') as img_file:
    for img_input in img_file:
        img_inputs.append(img_input[:-1])
#img_inputs = args.inputs
print(img_inputs[0])
print('len(img_inputs):',len(img_inputs))
assert len(img_inputs) > 0

print('Setting IMG as main TTree')
print('N IMG files:',len(img_inputs))
print('IMG file[0]:',img_inputs[0])
tree = ROOT.TChain("fevt/RHTree")
for fi in img_inputs:
    tree.Add(fi)
    #break
nEvts = tree.GetEntries()
print('N evts in IMG ntuple:',nEvts)

# Initialize output ntuple
if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)
file_out = ROOT.TFile.Open("%s/%s_mantuple.root"%(args.outdir, sample), "RECREATE")
tree_out = ROOT.TTree("mantuple", "mantuple")

#############################
# Initialize branches for the m_as as single floats
# [PyROOT boilerplate for single float branches]
ma0 = np.zeros(1, dtype='float32')
ma1 = np.zeros(1, dtype='float32')
tree_out.Branch('ma0', ma0, 'ma0/F')
tree_out.Branch('ma1', ma1, 'ma1/F')
#'''
mgg = np.zeros(1, dtype='float32')
tree_out.Branch('mgg', mgg, 'mgg/F')
sieie0 = np.zeros(1, dtype='float32')
sieie1 = np.zeros(1, dtype='float32')
tree_out.Branch('sieie0', sieie0, 'sieie0/F')
tree_out.Branch('sieie1', sieie1, 'sieie1/F')
chiso0 = np.zeros(1, dtype='float32')
chiso1 = np.zeros(1, dtype='float32')
tree_out.Branch('chiso0', chiso0, 'chiso0/F')
tree_out.Branch('chiso1', chiso1, 'chiso1/F')
hoe0 = np.zeros(1, dtype='float32')
hoe1 = np.zeros(1, dtype='float32')
tree_out.Branch('hoe0', hoe0, 'hoe0/F')
tree_out.Branch('hoe1', hoe1, 'hoe1/F')
pt0 = np.zeros(1, dtype='float32')
pt1 = np.zeros(1, dtype='float32')
tree_out.Branch('pt0', pt0, 'pt0/F')
tree_out.Branch('pt1', pt1, 'pt1/F')
bdt0 = np.zeros(1, dtype='float32')
bdt1 = np.zeros(1, dtype='float32')
tree_out.Branch('bdt0', bdt0, 'bdt0/F')
tree_out.Branch('bdt1', bdt1, 'bdt1/F')
r90 = np.zeros(1, dtype='float32')
r91 = np.zeros(1, dtype='float32')
tree_out.Branch('r90', r90, 'r90/F')
tree_out.Branch('r91', r91, 'r91/F')
iphi0 = np.zeros(1, dtype='float32')
iphi1 = np.zeros(1, dtype='float32')
tree_out.Branch('iphi0', iphi0, 'iphi0/F')
tree_out.Branch('iphi1', iphi1, 'iphi1/F')
phi0 = np.zeros(1, dtype='float32')
phi1 = np.zeros(1, dtype='float32')
tree_out.Branch('phi0', phi0, 'phi0/F')
tree_out.Branch('phi1', phi1, 'phi1/F')
ieta0 = np.zeros(1, dtype='float32')
ieta1 = np.zeros(1, dtype='float32')
tree_out.Branch('ieta0', ieta0, 'ieta0/F')
tree_out.Branch('ieta1', ieta1, 'ieta1/F')
eta0 = np.zeros(1, dtype='float32')
eta1 = np.zeros(1, dtype='float32')
tree_out.Branch('eta0', eta0, 'eta0/F')
tree_out.Branch('eta1', eta1, 'eta1/F')
energy0 = np.zeros(1, dtype='float32')
energy1 = np.zeros(1, dtype='float32')
tree_out.Branch('energy0', energy0, 'energy0/F')
tree_out.Branch('energy1', energy1, 'energy1/F')
energyCorr0 = np.zeros(1, dtype='float32')
energyCorr1 = np.zeros(1, dtype='float32')
tree_out.Branch('energyCorr0', energyCorr0, 'energyCorr0/F')
tree_out.Branch('energyCorr1', energyCorr1, 'energyCorr1/F')
phoiso0 = np.zeros(1, dtype='float32')
phoiso1 = np.zeros(1, dtype='float32')
tree_out.Branch('phoiso0', phoiso0, 'phoiso0/F')
tree_out.Branch('phoiso1', phoiso1, 'phoiso1/F')
tkiso0 = np.zeros(1, dtype='float32')
tkiso1 = np.zeros(1, dtype='float32')
tree_out.Branch('tkiso0', tkiso0, 'tkiso0/F')
tree_out.Branch('tkiso1', tkiso1, 'tkiso1/F')

nRecoPho = np.zeros(1, dtype='int32')
tree_out.Branch('nRecoPho', nRecoPho, 'nRecoPho/I')

run = np.zeros(1, dtype='int64')
tree_out.Branch('run', run, 'run/I')
lumi = np.zeros(1, dtype='int64')
tree_out.Branch('lumi', lumi, 'lumi/I')
evt = np.zeros(1, dtype='int64')
tree_out.Branch('evt', evt, 'evt/I')
pu = np.zeros(1, dtype='float32')
tree_out.Branch('pu', pu, 'pu/F')

ma3b30 = np.zeros(1, dtype='float32')
ma3b31 = np.zeros(1, dtype='float32')
tree_out.Branch('ma3b30', ma3b30, 'ma3b30/F')
tree_out.Branch('ma3b31', ma3b31, 'ma3b31/F')
mafc0 = np.zeros(1, dtype='float32')
mafc1 = np.zeros(1, dtype='float32')
tree_out.Branch('mafc0', mafc0, 'mafc0/F')
tree_out.Branch('mafc1', mafc1, 'mafc1/F')
#'''
branches = tree_out.GetListOfBranches()
branches = [br.GetName() for br in branches]
print(branches)

#############################
# Load mass regressor for inference
print('>> Loading regressor model')
def load_model_json(model, path):
    # Newer pytorch version save models as zip archive
    # that can't be loaded by older pytorch versions
    # Save model as json in newer pytorch, load json in older
    # then copy contents to an empty model state_dict
    # ref: https://stackoverflow.com/questions/64141188/how-to-load-checkpoints-across-different-versions-of-pytorch-1-3-1-and-1-6-x-u
    data_dict = OrderedDict()
    with open(path, 'r') as f:
        data_dict = json.load(f)
    own_state = model.state_dict()
    for k, v in data_dict.items():
        #print('Loading parameter:', k)
        if not k in own_state:
            print('Parameter', k, 'not found in own_state!!!')
        if type(v) == list or type(v) == int:
            v = torch.tensor(v)
        own_state[k].copy_(v)
    model.load_state_dict(own_state)
    print('.. model loaded from json.')
import torch_resnet_concat as networks
#epoch = 80
#neg_mass = 300
## AOD
#model_file = glob.glob('../EB_Pi0_massreg_Pt20To100_2017PU_ext3_tzfixed_wrapfix/MODELS/DoublePi0Pt20To100_m0To1600_pythia8_PU2017_genDR10_aodDR16_nPhoN_PhoNeg%dTo0_wgts_EBtzo25_AOD_m0o1.6_ResNet_blocks3_seedPos_FC128x0_MAEloss_lr0.0005_epochs80_n778k_run0/model_epoch%d_*.pkl'%(neg_mass, epoch))
#assert len(model_file) == 1
#model_file = model_file[0]
#'''
model_file = '/uscms/home/mba2012/nobackup/h2aa/CMSSW_10_5_0/src/h2aa/ntupleAnalysis/Models/model_epoch80_mae0.1906.pkl'
#model_file = 'models/resnet/precise300/model_epoch89_val_mae0.1838.json'
#model_file = 'models/resnet/precise300/model_epoch86_val_mae0.1872.json'
#model_file = 'models/resnet/precise300/model_epoch90_val_mae0.1889.json'
resnet = networks.ResNet(2, 3, [16, 32], 128, 0).cuda()
if '.json' in model_file:
    load_model_json(resnet, model_file)
else:
    resnet.load_state_dict(torch.load(model_file)['model']) # doesnt work with models from newer pytorch
    #resnet.load_state_dict(torch.jit.load(model_file)['model']) #doesnt work with models from newer pytorch
resnet.eval()
print(model_file)

# Load photon NN fcnet
#scaler = joblib.load('pi0_scaler.pkl')
#scaler = pickle.load(open('pi0_scaler_pickle.pkl', 'rb'))
# scaler wont load in lpc due to missing joblib python module(pickle doesnt work)
# use hardcoded values instead:
scale = np.array([2.83126097e-01, 2.70984447e-03, 1.31190771e+01, 1.38268298e+00,
                  2.54267599e+00, 3.46374138e+01, 1.96283358e-02, 2.31166868e-03,
                  6.91103320e-01, 2.85732555e-05, 6.58185195e-02])
mu = np.array([ 9.06606724e-01,  1.06437328e-02,  6.55752668e+00,  7.36942195e-01,
                2.01469760e+00,  6.70971420e+01,  2.56851846e-02,  9.73124830e-03,
                1.03052384e-03, -4.09735331e-07,  8.75230438e-01])
def scaler_transform(pho_id):
    return (np.array(pho_id)-mu)/scale
import torch_fcnet as fc_networks
#epoch = 27
#neg_mass = 400
#model_paths['phoid'] = '../EB_Pi0_massreg_Pt20To100_2017PU_ext3_tzfixed_PhoVars'
#model_file = glob.glob('%s/MODELS/DoublePi0Pt20To100_m0To1600_pythia8_PU2017_genDR10_recoDR16_nPhoN_PhoNeg%dTo0_wgts_phoId_rescale_m0o1.6_FC128x3_MAEloss_lr0.0005_epochs30_n778k/model_epoch%d_*.pkl'%(model_paths[algo], neg_mass, epoch)) # AOD
fcmodel_file = 'models/fcnet/model_epoch27_mae0.2957.pkl'
fcnet = fc_networks.FCNet(nodes=128, layers=3).cuda()
fcnet.load_state_dict(torch.load(fcmodel_file)['model'])
fcnet.eval()

eb_scale = 25.
m0_scale = 1.6
def transform_y(y):
    return y/m0_scale
def inv_transform(y):
    return y*m0_scale
def shapeEB(eb):
    return np.array(eb).reshape(1,170,360)

#############################
c, h = {}, {}
wd, ht = int(440*1), int(400*1)

'''
for b in [-0.92]:

    for s in [0.011]:

        for r in [1.e5]:

            for ptcut in [30., 40., 50., 1.e5]:

                #expt = 'bdtgt%.2f_sieielt%.3f_rcIsolt%.3f'%(b, s, r)
                expt = 'bdtgt%.2f_sieielt%.3f_rcIsolt%.3f_ptlt%.f'%(b, s, r, ptcut)
                #expt = 'bdtgt%.2f_sieielt%.3f'%(b, s)
                print('>> Expt:',expt)
'''

#k = 'iphi'
#c[k] = ROOT.TCanvas("c%s"%k,"c%s"%k,wd,ht)
#h[k] = ROOT.TH1F(k, k, 360//5, 0., 360)
#k = 'ieta'
#c[k] = ROOT.TCanvas("c%s"%k,"c%s"%k,wd,ht)
#h[k] = ROOT.TH1F(k, k, 170//5, 0., 170)

k = 'ma0'
#k += '-%s'%expt
##c[k] = ROOT.TCanvas("c%s"%k,"c%s"%k,wd,ht)
h[k] = ROOT.TH1F(k, k, 140, -0.2, 1.2)
#h[k] = ROOT.TH1F(k, k, 56, -0.2, 1.2)
##h[k] = ROOT.TH1F(k, k, 48, 0., 1.2)
k = 'ma1'
#k += '-%s'%expt
#c[k] = ROOT.TCanvas("c%s"%k,"c%s"%k,wd,ht)
h[k] = ROOT.TH1F(k, k, 140, -0.2, 1.2)
#h[k] = ROOT.TH1F(k, k, 56, -0.2, 1.2)
#h[k] = ROOT.TH1F(k, k, 48, 0., 1.2)

# Event range to process
iEvtStart = 0
iEvtEnd   = nEvts
#iEvtEnd   = 100#50000#10
print(">> Processing entries: [",iEvtStart,"->",iEvtEnd,")")

nWrite = 0
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    if iEvt%10e3==0: print(iEvt,'/',iEvtEnd-iEvtStart)
    evt_status = tree.GetEntry(iEvt)

    nPhoEvt = len(tree.SC_ieta)

    npho_roi = 0
    # Only keep events with photons within ieta image window
    for i in range(nPhoEvt):
        #if (tree.SC_ieta[i] < 15) or (tree.SC_ieta[i]+16 > 169): continue
        #if tree.pho_pT[i] > 100.: continue
        npho_roi += 1

    if npho_roi < 2: continue

    #if args.wgt_files is not None:
    #    rands = np.random.random((npho_roi, nPasses))

    #for i in range(npho_roi):
    #    h['ieta'].Fill(tree.SC_ieta[i])
    #    h['iphi'].Fill(tree.SC_iphi[i])

    #if tree.m0 > 110.: continue # GeV
    #if tree.pho_bdt[0] < 0.: continue
    #if tree.pho_bdt[1] < -0.98: continue
    #if tree.pho_sieie[1] > 0.012: continue
    #if tree.pho_chgIso[1]/tree.pho_pT[1] > 0.05: continue

    #if tree.pho_pT[1] > ptcut: continue
    #if tree.pho_bdt[1] < b: continue
    #if tree.pho_sieie[1] > s: continue
    #if tree.pho_chgIso[1]/tree.pho_pT[1] > r: continue
    #print('>> bdt: %f %f'%(tree.pho_bdt[0], tree.pho_bdt[1]))

    # Only keep events with barrel photons
    sc_cms = []
    ieta, iphi = [], []
    #X_EBt = np.array(tree.EB_energyT).reshape(1,170,360)
    #X_EBt = tree.EB_energyT
    #X_EBt = shapeEB(tree.EB_energyT)
    #print(type(tree.EB_energyT))
    #X_EBt = list(X_EBt)
    #del X_EBt
    #X_EBt = shapeEB(tree.EB_energyT)
    #X_EBz = shapeEB(tree.EB_energyZ)
    #'''
    X_EBt = np.array(tree.EB_energyT).reshape(1,170,360)
    X_EBz = np.array(tree.EB_energyZ).reshape(1,170,360)
    X_cms = np.concatenate([X_EBt, X_EBz], axis=0)
    #X_EBerr = np.array(tree.EB_energyErr).reshape(1,170,360)
    #print(X_cms.shape)
    for i in range(npho_roi):
        ieta.append([tree.SC_ieta[i]])
        iphi.append([tree.SC_iphi[i]])
        #print(tree.SC_ieta[pho_idx[i]], ieta[-1])
        #print(tree.SC_iphi[pho_idx[i]], iphi[-1])
        #sc_cms_ = crop_EBshower(X_cms, ieta[-1], iphi[-1])
        sc_cms_ = crop_EBshower_padded(X_cms, ieta[-1], iphi[-1])
        sc_cms.append(sc_cms_)
        #print(sc_cms_.shape)
    #print(np.array(sc_cms).shape)
    #print(np.array(iphi).shape)
    #print(np.array(ieta).shape)
    ma = inv_transform(resnet([torch.Tensor(sc_cms).cuda()/eb_scale,\
                               torch.Tensor(iphi).cuda()/360.,\
                               torch.Tensor(ieta).cuda()/170.
                              ])).tolist()
    #'''
    #print(np.array(ma).shape)
    #print(ma)
    #'''
    ma0_, ma1_ = ma[0][0], ma[1][0]
    #ma0_, ma0dn_, ma0up_, ma1_, ma1dn_, ma1up_ = np.array(ma)[:,0]
    #print(ma)
    #print('   .. ma0: %f, ma1: %f'%(ma0_, ma1_))
    ma0[0] = ma0_
    ma1[0] = ma1_
    #if tree.m0 > 110.: continue # GeV
    #if tree.pho_bdt[0] < 0.: continue
    #if tree.pho_bdt[1] < -0.98: continue
    #if tree.pho_sieie[1] > 0.012: continue
    #if tree.pho_chgIso[1]/tree.pho_pT[1] > 0.05: continue
    #'''
    mgg[0] = tree.m0
    sieie0[0] = tree.pho_sieie[0]
    sieie1[0] = tree.pho_sieie[1]
    chiso0[0] = tree.pho_chgIso[0]
    chiso1[0] = tree.pho_chgIso[1]
    hoe0[0] = tree.pho_HoE[0]
    hoe1[0] = tree.pho_HoE[1]
    pt0[0] = tree.pho_pT[0]
    pt1[0] = tree.pho_pT[1]
    bdt0[0] = tree.pho_bdt[0]
    bdt1[0] = tree.pho_bdt[1]

    r90[0] = tree.pho_r9[0]
    r91[0] = tree.pho_r9[1]
    iphi0[0] = tree.SC_iphi[0]
    iphi1[0] = tree.SC_iphi[1]
    phi0[0] = tree.pho_phi[0]
    phi1[0] = tree.pho_phi[1]
    ieta0[0] = tree.SC_ieta[0]
    ieta1[0] = tree.SC_ieta[1]
    eta0[0] = tree.pho_eta[0]
    eta1[0] = tree.pho_eta[1]
    energy0[0] = tree.pho_E[0]
    energy1[0] = tree.pho_E[1]
    energyCorr0[0] = -1 if '2018' in sample else tree.pho_ecalEPostCorr[0]
    energyCorr1[0] = -1 if '2018' in sample else tree.pho_ecalEPostCorr[1]

    phoiso0[0] = tree.pho_phoIso[0]
    phoiso1[0] = tree.pho_phoIso[1]
    tkiso0[0] = tree.pho_trkIso[0]
    tkiso1[0] = tree.pho_trkIso[1]

    nRecoPho[0] = tree.nRecoPho

    run[0] = tree.runId
    lumi[0] = tree.lumiId
    evt[0] = tree.eventId

    pu[0] = get_pu(tree.runId, tree.lumiId, inputPileupRange) if 'data' in sample else tree.evtPU
    #pu[0] = get_pu(tree.runId, tree.lumiId, inputPileupRange)
    #print('   .. run: %d, ls: %d, evt: %d, pu: %f'%(tree.runId, tree.lumiId, tree.eventId, pu[0]))

    pho_ids = []
    for i in range(npho_roi):
        pho_id_ = [
                tree.pho_r9[i]
                ,tree.pho_sieie[i]
                ,tree.pho_phoIso[i]
                ,tree.pho_chgIso[i]
                ,tree.pho_chgIsoWrongVtx[i]
                ,tree.pho_Eraw[i]
                ,tree.pho_phiWidth[i]
                ,tree.pho_etaWidth[i]
                ,tree.pho_scEta[i]
                ,tree.pho_sieip[i]
                ,tree.pho_s4[i]
        ]
        pho_ids.append(pho_id_)

    #print(np.array(pho_ids).shape)
    #print(pho_ids)
    pho_ids = scaler_transform(pho_ids)
    #print(np.array(pho_ids).shape)
    #print(pho_ids)

    mafc = inv_transform(fcnet(torch.Tensor(pho_ids).cuda())).tolist()
    mafc0_, mafc1_ = mafc[0][0], mafc[1][0]
    mafc0[0] = mafc0_
    mafc1[0] = mafc1_
    #print('   .. mafc0: %f, mafc1: %f'%(mafc0_, mafc1_))
    #'''
    ma3b30[0] = tree.pho_mass3b3[0]
    ma3b31[0] = tree.pho_mass3b3[1]

    tree_out.Fill()
    #h['ma0-%s'%expt].Fill(ma0_)
    #h['ma1-%s'%expt].Fill(ma1_)
    h['ma0'].Fill(ma0_)
    h['ma1'].Fill(ma1_)

    #break
    nWrite += 1

sw.Stop()
print(">> N events written: %d / %d"%(nWrite, iEvtEnd-iEvtStart))
print(">> Real time:",sw.RealTime()/60.,"minutes")
print(">> CPU time: ",sw.CpuTime() /60.,"minutes")

def set_hist(h, c, xtitle, ytitle, htitle):
    c.SetLeftMargin(0.16)
    c.SetRightMargin(0.15)
    c.SetBottomMargin(0.13)
    ROOT.gStyle.SetOptStat(0)

    h.GetXaxis().SetLabelSize(0.04)
    h.GetXaxis().SetLabelFont(62)
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleFont(62)
    h.GetXaxis().SetTitle(xtitle)

    h.GetYaxis().SetLabelSize(0.04)
    h.GetYaxis().SetLabelFont(62)
    h.GetYaxis().SetTitleOffset(1.5)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleFont(62)
    h.GetYaxis().SetTitle(ytitle)

    h.SetTitleSize(0.04)
    h.SetTitleFont(62)
    h.SetTitle(htitle)
    h.SetTitleOffset(1.2)

    return h, c

mass_line = 0.135

#######################################################
'''
#for shift in ['', 'dn', 'up']:
    k = 'ma0'+shift
    print(k, h.keys(), c.keys())
    assert k in h.keys()
    print(type(h[k]), type(c[k]))
    c[k].cd()
    h[k], c[k] = set_hist(h[k], c[k], "m_{a,pred} [GeV]", "N_{a}", "")
    h[k].SetLineColor(9)
    h[k].Draw("hist")
    k = 'ma1'+shift
    h[k].SetLineColor(2)
    h[k].Draw("hist SAME")

    ymax = 1.2*max(h['ma0'+shift].GetMaximum(), h['ma1'+shift].GetMaximum())
    h['ma0'+shift].GetYaxis().SetRangeUser(0., ymax)

    #l = ROOT.TLine(mass_line, 0., mass_line, ymax) # x0,y0, x1,y1
    #l.SetLineColor(14)
    #l.SetLineStyle(7)
    #l.Draw("same")
    hatch = ROOT.TGraph(2, array('d',[0.,0.]), array('d',[0.,ymax]));
    hatch.SetLineColor(14)
    hatch.SetLineWidth(2001)
    hatch.SetFillStyle(3004)
    hatch.SetFillColor(14)
    hatch.Draw("same")

    #c[k].SetGrid()
    c['ma0'+shift].Draw()
    c['ma0'+shift].Update()
    c['ma0'+shift].Print('Plots/h%s.eps'%k)
'''
'''
for im in ['ma0', 'ma1']:
    k = im
    c[k].cd()
    h[k], c[k] = set_hist(h[k], c[k], "m_{a,pred} [GeV]", "N_{a}", "")
    h[k].SetLineColor(1)
    h[k].Draw("hist")
    #k = im+'dn'
    #h[k].SetLineColor(2)
    #h[k].Draw("hist SAME")
    #k = im+'up'
    #h[k].SetLineColor(9)
    #h[k].Draw("hist SAME")

    ymax = 1.2*h[im].GetMaximum()
    h[im].GetYaxis().SetRangeUser(0., ymax)

    l = ROOT.TLine(mass_line, 0., mass_line, ymax) # x0,y0, x1,y1
    l.SetLineColor(14)
    l.SetLineStyle(7)
    l.Draw("same")
    hatch = ROOT.TGraph(2, array('d',[0.,0.]), array('d',[0.,ymax]));
    hatch.SetLineColor(14)
    hatch.SetLineWidth(2001)
    hatch.SetFillStyle(3004)
    hatch.SetFillColor(14)
    hatch.Draw("same")

    #c[k].SetGrid()
    c[im].Draw()
    c[im].Update()
    #c[im].Print('Plots/h%s_%s_nomdnup.eps'%(args.sample, im))
'''
'''
##############################
k = 'iphi'
c[k].cd()
h[k], c[k] = set_hist(h[k], c[k], "i#varphi", "N_{a}", "i#varphi")
h[k].SetLineColor(9)
h[k].Draw("hist")
c[k].Draw()
##############################
k = 'ieta'
c[k].cd()
h[k], c[k] = set_hist(h[k], c[k], "i#eta", "N_{a}", "i#eta")
h[k].SetLineColor(9)
h[k].Draw("hist")
c[k].Draw()
##############################

#blue: leading
#red: sub-leading
'''
file_out.Write()
file_out.Close()
