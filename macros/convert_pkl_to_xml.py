## to be used with
##
## python convert_pkl_to_xml.py models/2020Feb24/bdt_cmssw_mva_id_nnclean2/
## with the .pkl file in models/2020Feb24/bdt_cmssw_mva_id_nnclean2/
##
## At the end, the fixed.xml file has to be lint with
## xmllint --format filename.fixed.xml > filename.fixed.lint.xml
## this is the file to be used in RecoEgamma/ElectronIdentification

import os
from glob import glob
from argparse import ArgumentParser

import re
import xml.etree.cElementTree as ET

regex_float_pattern = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"

def build_tree(xgtree, base_xml_element, var_indices):
    parent_element_dict = {"0": base_xml_element}
    pos_dict = {"0": "s"}
    for line in xgtree.split("\n"):
        if not line:
            continue
        if ":leaf=" in line:

            result = re.match(r"(\t*)(\d+):leaf=({0})$".format(regex_float_pattern), line)
            if not result:
                print(line)
            depth = result.group(1).count("\t")
            inode = result.group(2)
            res = result.group(3)
            node_elementTree = ET.SubElement(
                parent_element_dict[inode],
                "Node",
                pos=str(pos_dict[inode]),
                depth=str(depth),
                NCoef="0",
                IVar="-1",
                Cut="0.0e+00",
                cType="1",
                res=str(res),
                rms="0.0e+00",
                purity="0.0e+00",
                nType="-99",
            )
        else:
            result = re.match(
                r"(\t*)([0-9]+):\[(?P<var>.+)<(?P<cut>{0})\]\syes=(?P<yes>\d+),no=(?P<no>\d+)".format(
                    regex_float_pattern
                ),
                line,
            )
            if not result:
                print(line)
            depth = result.group(1).count("\t")
            inode = result.group(2)
            var = result.group("var")
            cut = result.group("cut")
            lnode = result.group("yes")
            rnode = result.group("no")
            pos_dict[lnode] = "l"
            pos_dict[rnode] = "r"
            node_elementTree = ET.SubElement(
                parent_element_dict[inode],
                "Node",
                pos=str(pos_dict[inode]),
                depth=str(depth),
                NCoef="0",
                IVar=str(var_indices[var]),
                Cut=str(cut),
                cType="1",
                res="0.0e+00",
                rms="0.0e+00",
                purity="0.0e+00",
                nType="0",
            )
            parent_element_dict[lnode] = node_elementTree
            parent_element_dict[rnode] = node_elementTree

def convert_model(model, input_variables, output_xml):
    NTrees = len(model)
    var_list = input_variables
    var_indices = {}

    # <MethodSetup>                                                                                                    
    MethodSetup = ET.Element("MethodSetup", Method="BDT::BDT")

    # <Variables>                                                                                                                  
    Variables = ET.SubElement(MethodSetup, "Variables", NVar=str(len(var_list)))
    for ind, val in enumerate(var_list):
        name = val[0]
        var_type = val[1]
        var_indices[name] = ind
        Variable = ET.SubElement(
            Variables,
            "Variable",
            VarIndex=str(ind),
            Type=val[1],
            Expression=name,
            Label=name,
            Title=name,
            Unit="",
            Internal=name,
            Min="0.0e+00",
            Max="0.0e+00",
        )

    # <GeneralInfo>    
    GeneralInfo = ET.SubElement(MethodSetup, "GeneralInfo")
    Info_Creator = ET.SubElement(GeneralInfo, "Info", name="Creator", value="xgboost2TMVA")
    Info_AnalysisType = ET.SubElement(GeneralInfo, "Info", name="AnalysisType", value="Classification")

    # <Options>  
    Options = ET.SubElement(MethodSetup, "Options")
    Option_NodePurityLimit = ET.SubElement(Options, "Option", name="NodePurityLimit", modified="No").text = "5.00e-01"
    Option_BoostType = ET.SubElement(Options, "Option", name="BoostType", modified="Yes").text = "Grad"

    # <Weights>
    Weights = ET.SubElement(MethodSetup, "Weights", NTrees=str(NTrees), AnalysisType="1")

    for itree in range(NTrees):
        BinaryTree = ET.SubElement(Weights, "BinaryTree", type="DecisionTree", boostWeight="1.0e+00", itree=str(itree))
        build_tree(model[itree], BinaryTree, var_indices)

    tree = ET.ElementTree(MethodSetup)
    tree.write(output_xml)
    # format it with 'xmllint --format'   

        
parser = ArgumentParser()
parser.add_argument('model', default='$(pwd)/models/2020Feb24/bdt_cmssw_mva_id_nnclean2')
args = parser.parse_args()

print args.model

pkls = glob('%s/*.pkl' % args.model)
print pkls
if len(pkls) != 1 : raise RuntimeError('There must be one and only one pkl in the directory')

model = pkls[0]
what = 'cmssw_mva_id_nnclean2'

#this should be outsorced
from features import *
features, additional = get_features(what)
print 'Features'
print features
print 'Additional'
print additional

print 'Loading model'
from sklearn.externals import joblib
import xgboost as xgb
xml = model.replace('.pkl', '.xml')
model = joblib.load(model)

#from xgbo.xgboost2tmva import convert_model
from itertools import cycle

# xgb sklearn API assigns default names to the variables, use that to dump the XML
# then convert them to the proper name
print 'XML conversion'
xgb_feats = ['f%d' % i for i in range(len(features))]
print xgb_feats

convert_model(model._Booster.get_dump(), zip(xgb_feats, cycle('F')), xml)
xml_str = open(xml).read()
for idx, feat in reversed(list(enumerate(features))):
    xml_str = xml_str.replace('f%d' % idx, feat)

with open(xml.replace('.xml', '.fixed.xml'), 'w') as XML:
    XML.write(xml_str)
