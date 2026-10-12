"""Frozen gene sets (mouse symbols). Do not edit after pre-registration (notebook 00 hashes this file)."""
MITOPHAGY = ["Pink1", "Prkn", "Tbk1", "Optn", "Calcoco2", "Tax1bp1", "Bnip3", "Bnip3l", "Fundc1",
             "Fkbp8", "Bcl2l13", "Usp30", "Mul1", "Ambra1", "Vps13c", "Ulk1", "Rhot1", "Rhot2", "Mfn2"]
LYSOSOME = ["Atp6v0a1", "Atp6v0b", "Atp6v0c", "Atp6v0d1", "Atp6v1a", "Atp6v1b2", "Atp6v1c1", "Atp6v1d",
            "Atp6v1e1", "Atp6v1f", "Atp6v1h", "Atp6ap1", "Atp6ap2", "Ctsa", "Ctsb", "Ctsd", "Ctsl", "Ctsz",
            "Ctsc", "Gba", "Gla", "Gusb", "Hexa", "Hexb", "Man2b1", "Naglu", "Idua", "Sgsh", "Asah1",
            "Lipa", "Pla2g15", "Tpp1", "Ppt1", "Gns", "Galns", "Gaa", "Npc1", "Npc2", "Lamp1", "Lamp2",
            "Cd63", "Mcoln1", "Tmem175", "Slc38a9"]
AUTOPHAGY_CORE = ["Atg3", "Atg4b", "Atg5", "Atg7", "Atg9a", "Atg10", "Atg12", "Atg16l1", "Becn1",
                  "Map1lc3a", "Map1lc3b", "Gabarapl1", "Gabarapl2", "Sqstm1", "Ulk2", "Wipi2"]
PROTEASOME = ["Psma1", "Psma2", "Psma3", "Psma4", "Psma5", "Psma6", "Psma7", "Psmb1", "Psmb2", "Psmb3",
              "Psmb4", "Psmb5", "Psmb6", "Psmb7", "Psmc1", "Psmc2", "Psmc3", "Psmc4", "Psmc5", "Psmc6",
              "Psmd1", "Psmd2", "Psmd3", "Psmd4", "Psmd7", "Psmd8", "Psmd11", "Psmd12", "Psmd13", "Psmd14"]
RIBOSOME_PREFIX = ("Rpl", "Rps")
OXPHOS_PREFIX = ("Ndu", "Sdh", "Uqcr", "Cox4", "Cox5", "Cox6", "Cox7", "Atp5")
MICROGLIA = ["P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb", "Sall1", "Fcrls", "Aif1", "Trem2", "Tyrobp"]
NON_MICROGLIA = ["Snap25", "Rbfox3", "Slc17a7", "Gad1", "Plp1", "Mbp", "Mog", "Gfap", "Aqp4", "Pdgfra",
                 "Cldn5", "Pecam1", "Cd3e", "Cd79a", "Ms4a1"]
STATE_PANELS = {
    "homeostatic": ["P2ry12", "Tmem119", "Cx3cr1", "Selplg", "Siglech", "Hexb"],
    "DAM": ["Apoe", "Lpl", "Cst7", "Itgax", "Spp1", "Clec7a", "Axl", "Cd9", "Trem2", "Cd68"],
    "interferon": ["Ifit1", "Ifit2", "Ifit3", "Isg15", "Irf7", "Usp18", "Stat1", "Oasl2"],
    "chemokine_inflam": ["Ccl2", "Ccl3", "Ccl4", "Il1b", "Tnf", "Cxcl2", "Nfkbia"],
    "proliferating": ["Mki67", "Top2a", "Stmn1", "Birc5"],
    "MHC_II": ["Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"],
}
DISSOCIATION = ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Hspa1a", "Hspa1b", "Dusp1",
                "Zfp36", "Ier2", "Btg2", "Socs3", "Nr4a1"]
# Built-in mini regulon (TFEB-family 'CLEAR' network targets) used when decoupler is not installed.
CLEAR_TARGETS = ["Lamp1", "Ctsd", "Ctsb", "Ctsa", "Mcoln1", "Atp6v1h", "Atp6v0d1", "Atp6v1c1", "Gla", "Gba",
                 "Hexa", "Gns", "Sgsh", "Naglu", "Tpp1", "Ppt1", "Asah1", "Npc1", "Map1lc3b", "Sqstm1", "Becn1"]
