# Exact-solution references — provenance

Vendored 2026-09-04 from the co2-eos-cfd standalone (`suite/`), commit `ffc87e5beff22b22d8f7f0172d71712d1a019a4b`.
Generator: `suite/exact_riemann.py --case <B> --eos pr --model hem|frozen --N 800`
(`profiles/<case>.csv` are `--model frozen`; `exact_<B>_pr.csv` are HEM; `_frozen` are the bracket rows for B2/B9;
`_gerg` are the GERG-EOS variants used by convergence.py only).

Ground rule 19: these files have authority because they are computed independently of CAMR.
They are frozen — never regenerate them as part of a CAMR change. To regenerate after a
deliberate standalone change, re-run the generator, update this file, and re-baseline the gate.

Status in the standalone at vendoring time: `profiles/*.csv` are committed there;
the `exact_*.csv` files were UNTRACKED there (this copy is the only version-controlled one).

| file | bytes | mtime in standalone | md5 |
|---|---:|---|---|
| profiles/A1-Sod-strong.csv | 15359 | 2026-07-05 | 182e40dbf1b8feaab8aa9f2a49022bca |
| profiles/A2-Sod-weak.csv | 14994 | 2026-07-05 | a6af7e0d46b3e751341ee216a64e3776 |
| profiles/A3-Lax-like.csv | 15410 | 2026-07-05 | 92ffef85996d1899b9a731bcbd88499e |
| profiles/A4-Double-rare.csv | 15186 | 2026-07-05 | f196a6da434425b897febae49b61f7ce |
| profiles/A5-Two-shock.csv | 15597 | 2026-07-05 | a278f2a0f05df1c8ccb6fdb6dbc1b51b |
| profiles/A6-Near-vacuum.csv | 15094 | 2026-07-05 | 3f2fcaa37c8cad8f6354ce3a0b408f4b |
| profiles/B1-Comp-L-expand.csv | 15361 | 2026-07-05 | 57f0b98a41fb97354de3c6466d655f43 |
| profiles/B10-Cross-critical-hot.csv | 15203 | 2026-07-05 | d1ca32c8f6ef81ff355fffb9568f9671 |
| profiles/B2-Evap-wave.csv | 15061 | 2026-07-05 | 84bee27efd07f35e3c3f0dd30e19304d |
| profiles/B3-Sat-LV-contact.csv | 15286 | 2026-07-05 | f8f842f3f29b7cc85038a9101399024f |
| profiles/B4-Cross-critical.csv | 14958 | 2026-07-05 | d2a6244da001c0f4e33fd1fd3df5439b |
| profiles/B5-Both-2P.csv | 15293 | 2026-07-05 | 6b6d84d6715473705cfafd47d755e246 |
| profiles/B6-Sat-V-shock.csv | 15311 | 2026-07-05 | 73e35ab76cac679e7b058b509329e2ff |
| profiles/B7-Rupture-Sonic.csv | 15380 | 2026-07-05 | c5568933096084f5d340974d0261bb3d |
| profiles/B8-Wall-Reflection.csv | 15638 | 2026-07-05 | b54e046af42fdcf2cbf6e6a1f692b85b |
| profiles/B9-Deep-Expansion.csv | 15148 | 2026-07-05 | 0fa96032bc2613afd3b01fc391932318 |
| profiles/C1-Identity.csv | 15065 | 2026-07-05 | fb091b01665fbf577d9a9cd44fd1b037 |
| profiles/C2-Acoustic-limit.csv | 14670 | 2026-07-05 | ddb7f1c08bea5f06614937c981cb15fa |
| profiles/C3-Strong-shock-V.csv | 15293 | 2026-07-05 | b960d61d6d6e35caecc81433f34c7000 |
| exact_B10_pr.csv | 68014 | 2026-08-12 | a63e6eece47ea81d0b473713f3c390bd |
| exact_B11_pr.csv | 69082 | 2026-08-13 | c1fc0c23737364e90e01991b5d5c4392 |
| exact_B1_pr.csv | 68014 | 2026-08-08 | eece5e9f777313938fd09694271162a5 |
| exact_B2_pr.csv | 68014 | 2026-08-08 | eb3732a6a98a0b9096df9d19df1c7b7d |
| exact_B2_pr_frozen.csv | 68014 | 2026-08-15 | 4c8064147b61da4982c059f8c4d201b0 |
| exact_B3_pr.csv | 68014 | 2026-08-12 | 5aec404b68dd3efd777df9c597946ad5 |
| exact_B4_gerg.csv | 34014 | 2026-07-26 | b3a2c98b62b13c8300ea4b3b82b0fa6c |
| exact_B4_pr.csv | 68014 | 2026-07-26 | 63af1a0f870fad27a767c6a76214a88c |
| exact_B5_pr.csv | 68014 | 2026-08-12 | 076edf2809c379626b46d106d424ffc0 |
| exact_B6_pr.csv | 68014 | 2026-08-12 | 2dd22aa09c191ba2aac961865d191a3a |
| exact_B7_pr.csv | 68014 | 2026-08-12 | 9fefa0426fbd48f9b3e57e8793f078d6 |
| exact_B8_pr.csv | 68705 | 2026-08-12 | 33303641bcd3ffa9a0497d42b40314d2 |
| exact_B9_gerg.csv | 34014 | 2026-07-26 | a7e30cafa5ec4a777444c04c61eca099 |
| exact_B9_pr.csv | 68014 | 2026-07-26 | 5b8987dba4fd256efbff4cc110358007 |
| exact_B9_pr_frozen.csv | 68014 | 2026-08-15 | 82585b1299ac9fe5e6bc316b54f72e3c |
