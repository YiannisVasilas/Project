import cobra as cobra
import re
import json

def make_MNXRef_ns_dict(ns,reverse=False):
    ### PURPOSE
    # Creates a dictionary containing identifiers of namespace ns as keys and MNXRef identifiers as values
    ### EXAMPLE CALL
    #   name_dict = make_MNXRef_ns_dict('metacyc',reverse=False)
    ### INPUT
    #   ns: reference namespace. Current options are: metacyc, bigg, kegg, chebi, seed, upa, biopath, lipidmaps, hmdb, reactome, umbbd
    ### OPTIONAL INPUT
    #   reverse: If the dictionary should be swapped - MNXref ids as keys and ns ids as values
    ### OUTPUT
    #   name_dict - mapping dictionary
    
    #Input processing
    ns = ns.lower() #Remove caps
    
    #Create output structure
    name_dict = {}
    
    #makes library of reference compound id as key: mnx compound id as value
    with open('../Data/chem_xref.tsv') as f:
        for line in f:
            if line.startswith(ns+':'):
                elems = line.split('\t')
                elems[0] = re.sub('.*:','',elems[0])
                if reverse:
                    name_dict[elems[1]] = elems[0]
                else:
                    name_dict[elems[0]] = elems[1]
    
    #some manual additions / corrections:
    if ns=='bigg':
        name_dict['hbut4coa'] = 'MNXM2334'
        name_dict['mylcoa'] = 'MNXM1611'
        name_dict['mcnlcoa'] = 'MNXM2693'
        name_dict['mmylcoa'] = 'MNXM5993'
        name_dict['hbut4'] = 'MNXM514'
        name_dict['citmcoa'] = 'MNXM2358'
        name_dict['3hpcoa'] = 'MNXM1263'
        name_dict['2h3opp'] = 'MNXM475'
        name_dict['CPD1065'] = 'MNXM5823'
        name_dict['fdxo_4_2'] = 'MNXM145541'
        name_dict['fdxr_4_2'] = 'MNXM145542'
        name_dict['fdxr_2_2'] = 'MNXM169'
        name_dict['f6p'] = 'MNXM162235'
        name_dict['3oodcoa_c'] = 'MNXM513'
        name_dict['4cmcoa_c'] = 'MNXM264'
        name_dict['oxoarg__L'] = 'MNXM1037'
        name_dict['chpd_c'] = 'MNXM9829'
        name_dict['4h2oxov_c'] = 'MNXM1176'
        name_dict['5mthglu_c'] = ''
        name_dict['thglu_c'] = ''
        name_dict['quin_kt_c'] = ''
        name_dict['4iz5pp_c'] = ''
        name_dict['urocan_c'] = ''
        name_dict['slcys_c'] = ''
        name_dict['codscl2_c'] = ''
        name_dict['colipaA_c'] = ''
        name_dict['glc_bD_c'] = ''
        name_dict['glc_B_c'] = ''
        name_dict['bcarote_c'] = ''
        name_dict['atocophe_c'] = ''
        name_dict['btocophe_c'] = ''

    elif ns == 'seed':
        name_dict['cpd00027'] = 'MNXM41' 

    #Check if empty
    if len(name_dict)==0:
        print ('Unknown reference namespace. Empty dictionary returned')
    
    return name_dict

def Translate_to_MNXRef(model,ns):
    ### PURPOSE
    # Converts a model into the MNXref namespace
    ### EXAMPLE CALL
    #   Translate_to_MNXRef(model,'metacyc')
    ### INPUT
    #   model: CobraPy model structure
    #   ns: Current namespace in model
    ### OPTIONAL INPUT
    ### OUTPUT
    
    name_dict = make_MNXRef_ns_dict(ns)
    if len(name_dict)==0:
        print 'Could not map: unknown namespace'
        return
    
    print 'Mapping model with %d metabolites' % len(model.metabolites)
    print '%d metabolites require mapping' % len([m for m in model.metabolites if m.id[0:4]!='MNXM'])

    metlist = [m.id for m in model.metabolites]
    
    for mid in metlist:
        m = model.metabolites.get_by_id(mid)
        compound = m.id[:-2]
        compartment = m.id[-2:]
        if compound in name_dict:
            newid = name_dict[compound]+compartment
            #check for multiple metabolites ending up with the same name, remove them after replacing one of them with the other in its reactions:
            if newid not in model.metabolites:
                m.id = newid
                model.repair()
            else:
                replace_delete(model, m, model.metabolites.get_by_id(newid))
        
    #then replace compound name, formula and charge
    MNXRef_dict = MNXRef_met_props() #Gather metabolite properties
    for m in model.metabolites:
        mid = m.id[:-2]
        if mid in MNXRef_dict:
            m.name = MNXRef_dict[mid]['name']
            m.formula = MNXRef_dict[mid]['formula']
            m.charge = MNXRef_dict[mid]['charge']
    model.repair()
    
    print '%d unmapped metabolites remain' % (len([m for m in model.metabolites if m.id[0:4]!='MNXM']))
    
def MNXRef_met_props():
    ### PURPOSE
    # Creates a dictionary of dictionaries containing all MNXRef compound information
    ### EXAMPLE CALL
    #   MNXRef_dict = MNXRef_met_props()
    ### INPUT
    ### OPTIONAL INPUT
    ### OUTPUT
    #   MNXRef_dict: Dictionary of dictionaries. Primary key is compound id, secondary keys are name, formula, and charge.
    
    MNXRef_dict = {}
    
    with open('../Data/chem_prop.tsv') as f:
        for line in f:
            if line.startswith('MNXM'):
                elems = line.split('\t')
                if elems[3]!='' and elems[3]!='NA':
                    MNXRef_dict[elems[0]] = {'name':elems[1],
                                             'formula':elems[2],
                                             'charge':int(elems[3])}
    
    #some manual corrections:
    MNXRef_dict['MNXM318']['charge'] = 0
    MNXRef_dict['MNXM33']['charge'] = -2
    MNXRef_dict['MNXM89582']['formula'] = 'H2S' 
    MNXRef_dict['MNXM89582']['charge'] = 0

    return MNXRef_dict

def replace_delete(model, old_met, new_met):
    ### PURPOSE
    # Replaces all instances of one metabolite with another.
    ### EXAMPLE CALL
    #   replace_delete(coli_model, GLC, D-Glucose)
    ### INPUT
    #   model: CobraPy model instance
    #   old_met: metabolite to replace
    #   new_met: metabolite that takes over the reactions from old_met. Note: metabolite should already exist
    ### OPTIONAL INPUT
    ### OUTPUT
    #   Updated model.
    
    reactions = [r for r in old_met.reactions]
    
    for r in reactions:
        replace_stochiometry = r.metabolites[old_met]
        r.add_metabolites({new_met: replace_stochiometry})
    #finally remove the replaced metabolite from the model to get it out of the reaction
    old_met.remove_from_model()
    model.repair()
    return

# def check_Keq():
#     REACTION_FNAME = '../Data/pathwayreactions.txt'
#     ccCalc = cc.component_contribution_trainer.ComponentContribution.init()
#     reaction_strings = open(REACTION_FNAME, 'r').readlines()
#     model = cc.kegg_model.KeggModel.from_formulas(reaction_strings)

#     model.add_thermo(ccCalc)
#     dG0_prime, dG0_std,  sqrt_sigma = model.get_transformed_dG0(7.5, 0.1, 298.15)
#     Keq = []
#     for i in dG0_prime:
#         Keq.append(math.exp(-1*i/(8.314E-3*298.15)))
#     return Keq

def mnx_to_keggformat(model):
    for m in model.metabolites:
        if m.id[-2:] == '_c' or m.id[-2:] == '_m':
            m.id = m.id[:-2]
    model.repair()
    mnxkegg = mnxkeggdict()
    for m in model.metabolites:
        if m.id in mnxkegg:
            mid_old = m.id
            if mnxkegg[m.id] in model.metabolites:
                replace_delete(model,m,model.metabolites.get_by_id(mnxkegg[m.id]))
            else:
                m.id = mnxkegg[m.id]
                model.repair()
    return

def kegg_reactions():
    #create dict RID: reaction
    f = open(equilibrator_loc+'data/kegg_reactions.json', 'r')
    d = json.load(f)
    result = {}
    for elem in d:
        result[elem['RID']] = elem['reaction']
    #get strings instead of unicode:
    result2 = {}
    for elem in result:
        key = str(elem)
        value = []
        result_value = result[elem]
        for m in result_value:
            m[1] = str(m[1])
            value.append(m)
        result2[key] = value
    return result2

def mnxkeggdict():
    #makes mnx - > kegg dictionary (for delta G determination)
    xref = open("../Data/chem_xref.tsv",  'r')
    mnxkegg = {}
    while True:
        line = xref.readline()
        if line[0:6] == 'kegg:C':
            tablocs = []
            for i in range(len(line)):
                if line[i] == '\t' or line[i] == '\n':
                    tablocs.append(i)
            if line[tablocs[0]+1:tablocs[1]] not in mnxkegg:
                mnxkegg[line[tablocs[0]+1:tablocs[1]]] = line[5:tablocs[0]]
        if line[0:10] == 'lipidmaps:':
            break
    #some manual additions:
    mnxkegg.update({'MNXM84':'C00014', 'MNXM191': 'C00390',  'MNXM232': 'C00399', 'MNXM145542': 'C00138', 'MNXM145541':'C00139',  'MNXM169': 'C00138', 'MNXM178':'C00139', 
    'MNXM145809':'C02869', 'MNXM6271':'C02745', 'MNXM148165':'C06021', 'MNXM148349':'C06020', 'MNXM558':'C00390', 'MNXM2178': 'C00399', 'MNXM509': 'C00828',  'MNXM223':'C05819', 
    'MNXM116':'C00117', 'MNXM15900':'C00117', 'MNXM162':'C00311','MNXM1783':'C00399',  'MNXM5319':'C00390', 'MNXM114366':'C14180','MNXM89629':'C00085','MNXM2378':'C00118','MNXM145523':'C00051','MNXM91032':'C18026',}) #Why 2? - '2 C00138', 'MNXM145541': '2 C00139',  'MNXM169': '2 C00138', 'MNXM178':'2 C00139'
    xref.close()
    return mnxkegg

def getkeggformationenergies(filename):
    # load formation energies from the JSON file
    COMPOUND_DICT = {}
    for cd in json.load(open(filename, 'rb')):
        kegg_id = cd.get('CID', 'unknown')
        COMPOUND_DICT[kegg_id] = cd

    return COMPOUND_DICT

def check_samereaction(reactiona,  reactionb):
    #checks if two reactions are exactly the same
    alist = []
    blist = []
    clist = []
    for elem in reactiona.metabolites:
        alist.append((elem.id,  reactiona.metabolites[elem]))
    for elem in reactionb.metabolites:
        blist.append((elem.id, reactionb.metabolites[elem]))
    #also check with exactly reverse reaction:
    for elem in reactionb.metabolites:
        clist.append((elem.id,  reactionb.metabolites[elem]*-1))
    same = True
    reverse = False
    asame1 = True
    asame2 = True
    bsame = True
    csame = True
    for elem in alist:
        if elem not in blist:
            bsame = False
        if elem not in clist:
            csame = False
    for elem in blist:
        if elem not in alist:
            asame1 = False
    for elem in clist:
        if elem not in alist:
            asame2 = False
    if asame1 == False and asame2 == False:
        same = False
    if bsame == False and csame == False:
        same = False
    if csame == True and asame2 == True:
        reverse = True
    return same#,  reverse

def silence_oxygen_sensitive(model):
    knockoutlist = ['POR5', 'POR_2', 'PFL', 'CODH_ACS', 'FDH7', 'RNF', 'OOR2', 'MNXR83963', 'MNXR4175', 
    'MNXR82027', 'MNXR82814', 'MNXR83529', 'MNXR83965', 'MNXR35431', 'MNXR54890', 'MNXR85543']
    for r in knockoutlist:
        if r in model.reactions:
            model.reactions.get_by_id(r).lower_bound = 0
            model.reactions.get_by_id(r).upper_bound = 0
    return
