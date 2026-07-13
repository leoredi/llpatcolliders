#!/usr/bin/env wolframscript
(* Decode the pinned SensCalc 2501.04525 MX files into reviewable formats. *)

ClearAll[
  die, render, plainString, canonicalDecayName, normalizedProducts, sha256,
  checkedImport, evaluateAtMass, writeJSON
];

die[message_] := (
  WriteString[$Output, "error: " <> ToString[message, InputForm] <> "\n"];
  Exit[1]
);

render[expression_] := ToString[expression, InputForm, PageWidth -> Infinity];
plainString[expression_] := If[StringQ[expression], expression, render[expression]];
canonicalDecayName[expression_] := Switch[plainString[expression],
  "2e", "ee",
  "2\[Mu]", "mumu",
  "2\[Tau]", "tautau",
  "2\[Gamma]", "gammagamma",
  "Non-hadronic, total", "nonhadronic_total",
  "Hadronic, total", "hadronic_total",
  "Total", "total",
  _, Null
];
normalizedProducts[products_] := plainString /@
  DeleteCases[Flatten[{products}], Null | "Null"];
sha256[path_] := IntegerString[FileHash[path, "SHA256"], 16, 64];

checkedImport[path_] := Module[{value},
  value = Quiet[Check[Import[path, "MX"], $Failed]];
  If[value === $Failed, die["failed to import MX file " <> path]];
  value
];

evaluateAtMass[expression_, mass_] := Module[{value},
  value = expression /. ((symbol_Symbol /;
       MemberQ[{"ma", "mLLP", "mN"}, SymbolName[Unevaluated[symbol]]]) :> mass);
  If[MatchQ[value, _Function | _InterpolatingFunction], value = value[mass]];
  value = Quiet[Check[N[value], Indeterminate]];
  If[!NumericQ[value] || !TrueQ[Im[N[value]] == 0],
    die["non-numeric branching ratio at mass " <> ToString[mass] <>
      ": " <> render[expression]]
  ];
  Re[N[value]]
];

writeJSON[path_, value_] := Module[{result},
  result = Quiet[Check[Export[path, value, "RawJSON"], $Failed]];
  If[result === $Failed, die["failed to write " <> path]]
];

If[Length[$ScriptCommandLine] != 3,
  die["usage: export_senscalc_2501.wls SENSCALC_ROOT OUTPUT_DIR"]
];

sensCalcRoot = ExpandFileName[$ScriptCommandLine[[2]]];
outputDir = ExpandFileName[$ScriptCommandLine[[3]]];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];

sourceFiles = <|
  "widths" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "decay widths",
    "Widths-model-ALP-fermion-scale-1000.-GeV-2501.04525.m"
  }],
  "branching_ratios" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "decay widths",
    "Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-2501.04525.m"
  }],
  "matrix_elements" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "decay widths",
    "Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-2501.04525.m"
  }],
  "acceptance_notebook" -> FileNameJoin[{
    "codes", "Acceptances", "ALP-fermion.nb"
  }]
|>;
expectedHashes = <|
  "widths" -> "d12fb78d28edff0aa081c9fb66d829b42b4ec71202684019d7b9047ecb40b869",
  "branching_ratios" -> "34a09eed87d081bfffe79b464094741454022f79478c5b28bc236bc361049286",
  "matrix_elements" -> "f959c2257fa349e5af6966795db8cbf0da2e3ff8d097b5fccc3316d270ff17ed",
  "acceptance_notebook" -> "8d205b65da456fb1a8fac8d1803eaa73ae839128a60d0624eeea17cacee5d4b4"
|>;
sourcePaths = AssociationMap[FileNameJoin[{sensCalcRoot, #}] &, sourceFiles];

KeyValueMap[
  Function[{key, path},
    If[!FileExistsQ[path], die["missing source file " <> path]];
    If[sha256[path] =!= expectedHashes[key], die["SHA-256 mismatch for " <> path]]
  ],
  sourcePaths
];

widthRaw = checkedImport[sourcePaths["widths"]];
If[!ListQ[widthRaw] || Length[widthRaw] < 4 || !ListQ[widthRaw[[4]]],
  die["unexpected widths MX structure"]
];
widthPayload = widthRaw[[4]];
If[Length[widthPayload] < 4, die["widths payload has no numeric rows"]];
widthSourceHeaders = widthPayload[[1]];
widthCanonicalNames = canonicalDecayName /@ widthSourceHeaders;
widthSourceMetadata = Take[widthPayload, 3];
widthDataRaw = N[Drop[widthPayload, 3]];
If[!MatrixQ[widthDataRaw, NumericQ], die["widths payload is not a numeric matrix"]];
If[Length[widthSourceHeaders] =!= Length[First[widthDataRaw]],
  die["width header and data-column counts differ"]
];

widthColumnIds = Join[
  {"mass_GeV"},
  Table["width_coefficient_" <> IntegerString[index, 10, 3],
    {index, Length[First[widthDataRaw]] - 1}]
];
requiredWidthNames = {
  "ee", "mumu", "tautau", "gammagamma", "nonhadronic_total",
  "hadronic_total", "total"
};
If[!And @@ (MemberQ[widthCanonicalNames, #] & /@ requiredWidthNames),
  die["required canonical width columns are missing"]
];
widthDataBNT = (Prepend[Rest[#]/4., First[#]] &) /@ widthDataRaw;
Export[
  FileNameJoin[{outputDir, "widths_raw_senscalc.csv"}],
  Prepend[widthDataRaw, widthColumnIds],
  "CSV"
];
Export[
  FileNameJoin[{outputDir, "widths_bnt.csv"}],
  Prepend[widthDataBNT, widthColumnIds],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "widths_metadata.json"}],
  <|
    "column_ids" -> widthColumnIds,
    "source_headers" -> (plainString /@ widthSourceHeaders),
    "source_headers_input_form" -> (render /@ widthSourceHeaders),
    "columns" -> Table[
      <|
        "source_index" -> index,
        "id" -> widthColumnIds[[index]],
        "source_header" -> plainString[widthSourceHeaders[[index]]],
        "canonical_name" -> widthCanonicalNames[[index]]
      |>,
      {index, Length[widthColumnIds]}
    ],
    "source_metadata_rows_input_form" -> (render /@ widthSourceMetadata),
    "raw_coefficient_definition" -> "Gamma = (gY/(2 vH))^2 * coefficient",
    "bnt_coefficient_definition" -> "Gamma = (1/f_BNT)^2 * coefficient",
    "conversion_raw_to_bnt" -> 0.25
  |>
];

branchingRaw = checkedImport[sourcePaths["branching_ratios"]];
If[!ListQ[branchingRaw] || Length[branchingRaw] < 3,
  die["unexpected branching-ratio MX structure"]
];
branchingProcesses = branchingRaw[[1]];
branchingProducts = branchingRaw[[2]];
branchingExpressions = branchingRaw[[3]];
If[!SameQ[
    Length[branchingProcesses], Length[branchingProducts],
    Length[branchingExpressions]
  ],
  die["branching process, product, and expression counts differ"]
];

channelIds = Table[
  "channel_" <> IntegerString[index, 10, 3],
  {index, Length[branchingProcesses]}
];
channelCanonicalNames = canonicalDecayName /@ branchingProcesses;
If[!And @@ (MemberQ[channelCanonicalNames, #] & /@
    {"ee", "mumu", "tautau", "gammagamma"}),
  die["required canonical branching channels are missing"]
];
massGrid = DeleteDuplicates[widthDataRaw[[All, 1]]];
branchingData = Table[
  Prepend[evaluateAtMass[#, mass] & /@ branchingExpressions, mass],
  {mass, massGrid}
];
Export[
  FileNameJoin[{outputDir, "branching_ratios.csv"}],
  Prepend[branchingData, Join[{"mass_GeV"}, channelIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "decay_channels.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> channelIds[[index]],
      "canonical_name" -> channelCanonicalNames[[index]],
      "process_name" -> plainString[branchingProcesses[[index]]],
      "decay_products" -> normalizedProducts[branchingProducts[[index]]],
      "process_input_form" -> render[branchingProcesses[[index]]],
      "decay_products_input_form" -> render[branchingProducts[[index]]],
      "branching_expression_input_form" -> render[branchingExpressions[[index]]]
    |>,
    {index, Length[channelIds]}
  ]
];

matrixRaw = checkedImport[sourcePaths["matrix_elements"]];
If[!ListQ[matrixRaw] || Length[matrixRaw] < 3,
  die["unexpected matrix-element MX structure"]
];
matrixProcesses = matrixRaw[[1]];
matrixMetadata = matrixRaw[[2]];
matrixExpressions = matrixRaw[[3]];
If[!SameQ[
    Length[matrixProcesses], Length[matrixMetadata],
    Length[matrixExpressions]
  ],
  die["matrix-element process, metadata, and expression counts differ"]
];
matrixIds = Table[
  "matrix_element_" <> IntegerString[index, 10, 3],
  {index, Length[matrixProcesses]}
];
writeJSON[
  FileNameJoin[{outputDir, "matrix_elements.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> matrixIds[[index]],
      "canonical_name" -> canonicalDecayName[matrixProcesses[[index]]],
      "process_name" -> plainString[matrixProcesses[[index]]],
      "decay_products" -> normalizedProducts[matrixMetadata[[index]]],
      "process_input_form" -> render[matrixProcesses[[index]]],
      "metadata_input_form" -> render[matrixMetadata[[index]]],
      "expression_input_form" -> render[matrixExpressions[[index]]],
      "expression_c_form" -> ToString[
        CForm[matrixExpressions[[index]]], PageWidth -> Infinity
      ]
    |>,
    {index, Length[matrixIds]}
  ]
];

outputFiles = {
  "widths_raw_senscalc.csv", "widths_bnt.csv", "widths_metadata.json",
  "branching_ratios.csv", "decay_channels.json", "matrix_elements.json"
};
writeJSON[
  FileNameJoin[{outputDir, "EXPORT_MANIFEST.json"}],
  <|
    "senscalc_repository" -> "https://github.com/maksymovchynnikov/SensCalc",
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
  " SensCalc data files to " <> outputDir <> "\n"];
