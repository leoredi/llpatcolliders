#!/usr/bin/env wolframscript
(* Decode the pinned SensCalc 2501.04525 LHC production inputs. *)

ClearAll[
  die, render, plainString, sha256, checkedImport, numericAtMass, writeJSON
];

die[message_] := (
  WriteString[$Output, "error: " <> ToString[message, InputForm] <> "\n"];
  Exit[1]
);

render[expression_] := ToString[expression, InputForm, PageWidth -> Infinity];
plainString[expression_] := If[StringQ[expression], expression, render[expression]];
sha256[path_] := IntegerString[FileHash[path, "SHA256"], 16, 64];

checkedImport[path_] := Module[{value},
  value = Quiet[Check[Import[path, "MX"], $Failed]];
  If[value === $Failed, die["failed to import MX file " <> path]];
  value
];

numericAtMass[expression_, mass_] := Module[{value},
  value = expression /. ((symbol_Symbol /;
       MemberQ[{"ma", "mLLP"}, SymbolName[Unevaluated[symbol]]]) :> mass);
  If[MatchQ[value, _Function | _InterpolatingFunction], value = value[mass]];
  value = Quiet[Check[N[value], Indeterminate]];
  If[!NumericQ[value] || !TrueQ[Abs[Im[N[value]]] < 10^-10],
    die["non-numeric production input at mass " <> ToString[mass] <>
      ": " <> render[expression]]
  ];
  Re[N[value]]
];

writeJSON[path_, value_] := Module[{result},
  result = Quiet[Check[Export[path, value, "RawJSON"], $Failed]];
  If[result === $Failed, die["failed to write " <> path]]
];

If[Length[$ScriptCommandLine] != 3,
  die["usage: export_senscalc_2501_production.wls SENSCALC_ROOT OUTPUT_DIR"]
];

sensCalcRoot = ExpandFileName[$ScriptCommandLine[[2]]];
outputDir = ExpandFileName[$ScriptCommandLine[[3]]];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];

sourceFiles = <|
  "experiments" -> FileNameJoin[{"codes", "experiments.nb"}],
  "coefficients" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "Coefficients-ALP-fermion-Lambda=1.-TeV.mx"
  }],
  "fragmentation_probability" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "ProductionProbability-Fragmentation-ALP-fermion.mx"
  }],
  "light_meson_decays" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "BrRatios-Msquared-LightMesonDecays-ALP-fermion-Lambda=1.-TeV.mx"
  }],
  "bremsstrahlung_lhc" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "ProdProb_ALP-fermion_Bremsstrahlung-AP_LHC.m"
  }],
  "drell_yan_lhc" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "sigmaDrellYan_ALP-fermion_LHC.txt"
  }],
  "fragmentation_lhc_grid" -> FileNameJoin[{
    "spectra", "New physics particles spectra", "ALP-fermion", "Pregenerated",
    "DoubleDistr_ALP-fermion_Fragmentation_LHC.m"
  }],
  "eta_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_Eta.txt"
  }],
  "eta_prime_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_EtaPr.txt"
  }],
  "omega_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_Omega.txt"
  }],
  "rho_charged_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_RhoCh.txt"
  }],
  "ks_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_KS.txt"
  }]
|>;

expectedHashes = <|
  "experiments" -> "3e60809105f49d8c5274de9a352cd0063bac9cb8feba31188e5e49ca6e468bdf",
  "coefficients" -> "e6885569ee8d33fb5310723ba2f020bf8b4e3506c2b5ad258b5230f9aa981f43",
  "fragmentation_probability" -> "4804205cfb14c0d40322fcabb7ec65843023ee44e91f6918603bca0e8835e9eb",
  "light_meson_decays" -> "c7468dcf8863a8ac182c2720d34d12b2dd12693700aec032460cd5d678985f71",
  "bremsstrahlung_lhc" -> "435022c25a2bb908ffd5c4b12a5645a0d5aa7f6914c49aab5e33fe0b8b60e3d2",
  "drell_yan_lhc" -> "12e013f6df9a65c07514e75d8f1c47d53b53ebc72b21624ab13eb70afaed185d",
  "fragmentation_lhc_grid" -> "43116651d526358c446af54a16d43ca9b2233fb42ff7ab4bae6946fe21477a79",
  "eta_lhc_spectrum" -> "a7e2323c721ce2f6b4b34571b9eaa956ce2cbb5ec8d5c437a2482cafdbfbf39c",
  "eta_prime_lhc_spectrum" -> "dfe394e4901f9d3ddf1e2e3881abe2f8fe550a1cbd0bc7a759f87b79fdf8882d",
  "omega_lhc_spectrum" -> "86f8a67260340de8091936a602eed44f33c5d78200afbaca5b5943fb5849eb02",
  "rho_charged_lhc_spectrum" -> "b2265c828c8b98762e6b7d99973ed6a97f5f66b01c2cab87708236f9ba4738ff",
  "ks_lhc_spectrum" -> "53aafc2e551fed9d01698ba37c7f230e6c5e66c30cc41e73976cd06dd21784ef"
|>;

sourcePaths = AssociationMap[FileNameJoin[{sensCalcRoot, #}] &, sourceFiles];
KeyValueMap[
  Function[{key, path},
    If[!FileExistsQ[path], die["missing source file " <> path]];
    If[sha256[path] =!= expectedHashes[key], die["SHA-256 mismatch for " <> path]]
  ],
  sourcePaths
];

auditMassGrid = N[Join[Range[0.02, 0.20, 0.01], Range[0.22, 1.00, 0.02]]];

coefficientsRaw = checkedImport[sourcePaths["coefficients"]];
writeJSON[
  FileNameJoin[{outputDir, "coefficients_metadata.json"}],
  <|
    "dimensions" -> Dimensions[coefficientsRaw],
    "input_form" -> render[coefficientsRaw]
  |>
];

fragmentationRaw = checkedImport[sourcePaths["fragmentation_probability"]];
If[!ListQ[fragmentationRaw] ||
   !And @@ ((ListQ[#] && Length[#] >= 2 &) /@ fragmentationRaw),
  die["unexpected fragmentation-probability MX structure"]
];
fragmentationFacilities = plainString /@ fragmentationRaw[[All, 1]];
fragmentationExpressions = fragmentationRaw[[All, 2]];
fragmentationIds = Table[
  "fragmentation_" <> IntegerString[index, 10, 3],
  {index, Length[fragmentationFacilities]}
];
fragmentationData = Table[
  Prepend[numericAtMass[#, mass] & /@ fragmentationExpressions, mass],
  {mass, auditMassGrid}
];
Export[
  FileNameJoin[{outputDir, "fragmentation_probability.csv"}],
  Prepend[fragmentationData, Join[{"mass_GeV"}, fragmentationIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "fragmentation_probability_metadata.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> fragmentationIds[[index]],
      "facility" -> fragmentationFacilities[[index]],
      "expression_input_form" -> render[fragmentationExpressions[[index]]]
    |>,
    {index, Length[fragmentationIds]}
  ]
];

lightRaw = checkedImport[sourcePaths["light_meson_decays"]];
If[!ListQ[lightRaw] ||
   !And @@ ((ListQ[#] && Length[#] >= 2 &) /@ lightRaw),
  die["unexpected light-meson MX structure"]
];
lightProcesses = plainString /@ lightRaw[[All, 1]];
lightBranchingExpressions = lightRaw[[All, 2]];
lightMatrixExpressions = (If[Length[#] >= 3, #[[3]], 0.] &) /@ lightRaw;
lightIds = Table[
  "light_meson_" <> IntegerString[index, 10, 3],
  {index, Length[lightProcesses]}
];
lightData = Table[
  Prepend[numericAtMass[#, mass] & /@ lightBranchingExpressions, mass],
  {mass, auditMassGrid}
];
Export[
  FileNameJoin[{outputDir, "light_meson_branching_ratios.csv"}],
  Prepend[lightData, Join[{"mass_GeV"}, lightIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "light_meson_channels.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> lightIds[[index]],
      "process_name" -> lightProcesses[[index]],
      "branching_expression_input_form" ->
        render[lightBranchingExpressions[[index]]]
    |>,
    {index, Length[lightIds]}
  ]
];
writeJSON[
  FileNameJoin[{outputDir, "light_meson_matrix_elements.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> lightIds[[index]],
      "process_name" -> lightProcesses[[index]],
      "expression_input_form" -> render[lightMatrixExpressions[[index]]],
      "expression_c_form" -> ToString[
        CForm[lightMatrixExpressions[[index]]], PageWidth -> Infinity
      ]
    |>,
    {index, Length[lightIds]}
  ]
];

fragmentationGrid = checkedImport[sourcePaths["fragmentation_lhc_grid"]];
If[!MatrixQ[fragmentationGrid, NumericQ] || Length[First[fragmentationGrid]] =!= 4,
  die["unexpected pregenerated fragmentation-grid structure"]
];
Export[
  FileNameJoin[{outputDir, "fragmentation_lhc_grid.csv"}],
  Prepend[N[fragmentationGrid], {"mass_GeV", "theta_rad", "energy_GeV", "density"}],
  "CSV"
];

drellYan = Quiet[Check[Import[sourcePaths["drell_yan_lhc"], "Table"], $Failed]];
If[drellYan === $Failed || !MatrixQ[drellYan, NumericQ] ||
   Length[First[drellYan]] =!= 4,
  die["unexpected Drell-Yan table structure"]
];
Export[
  FileNameJoin[{outputDir, "drell_yan_lhc.csv"}],
  Prepend[N[drellYan], {
    "mass_GeV", "sigma_central_pb_GeV2", "sigma_lower_pb_GeV2",
    "sigma_upper_pb_GeV2"
  }],
  "CSV"
];

bremsstrahlungRaw = checkedImport[sourcePaths["bremsstrahlung_lhc"]];
writeJSON[
  FileNameJoin[{outputDir, "bremsstrahlung_lhc_metadata.json"}],
  <|
    "dimensions" -> Dimensions[bremsstrahlungRaw],
    "input_form" -> render[bremsstrahlungRaw],
    "status" -> "audit-only; no transverse GRENDEL acceptance is implied"
  |>
];

outputFiles = {
  "coefficients_metadata.json", "fragmentation_probability.csv",
  "fragmentation_probability_metadata.json",
  "light_meson_branching_ratios.csv", "light_meson_channels.json",
  "light_meson_matrix_elements.json", "fragmentation_lhc_grid.csv",
  "drell_yan_lhc.csv", "bremsstrahlung_lhc_metadata.json"
};
writeJSON[
  FileNameJoin[{outputDir, "PRODUCTION_EXPORT_MANIFEST.json"}],
  <|
    "senscalc_tag" -> "v.1.3.3",
    "senscalc_commit" -> "0bca050633aae16e148d47f21840fa07ff4b8724",
    "phenomenology" -> "arXiv:2501.04525",
    "wolfram_version" -> $Version,
    "wolfram_system_id" -> $SystemID,
    "exported_utc" -> DateString[TimeZoneConvert[Now, 0], "ISODateTime"],
    "source_sha256" -> AssociationMap[sha256, sourcePaths],
    "output_sha256" -> AssociationMap[
      sha256[FileNameJoin[{outputDir, #}]] &,
      outputFiles
    ]
  |>
];

WriteString[$Output, "Exported " <> ToString[Length[outputFiles]] <>
  " SensCalc production data files to " <> outputDir <> "\n"];
