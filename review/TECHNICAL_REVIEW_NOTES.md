1|1|1|# rhoprime review notes
2|2|2|
3|3|3|## 2026-06-04 scope
4|4|4|User goal: act as a scientific-code foreman: understand the project, audit numerical/physics robustness, verify result reliability, and improve event-generator efficiency, especially unweighted/no-weight event generation.
5|5|5|
6|6|6|## Project summary
7|7|7|This repository models UPC production of rho'(1450) -> rho0 sigma -> pi+ pi- pi+ pi- in heavy-ion collisions. The flow is:
8|8|8|
9|9|9|1. `generate_grid_opt` computes d2sigma/dMdy and px/py spin-density histograms from a b-integral UPC model.
10|10|10|2. `generate_bose` samples M,y,px,py from those histograms and generates sequential four-pion decays with optional Bose symmetrization.
11|11|11|3. `analyze_bose` / `analyze_nonbose` reconstruct mass and azimuthal/spin-alignment observables.
12|12|12|
13|13|13|Main code evidence:
14|14|14|- `src/UpcRhoPrimeModelOpt.cxx`: optimized UPC grid and b integration.
15|15|15|- `src/EventGeneratorBose.cxx`: event generation, Bose amplitude, weights.
16|16|16|- `src/main_generate_bose.cxx`: CLI for weighted/unweighted event generation.
17|17|17|
18|18|18|## Baseline build and smoke tests
19|19|19|Build passed:
20|20|20|
21|21|21|```bash
22|22|22|cmake --build build -j$(nproc)
23|23|23|```
24|24|24|
25|25|25|Weighted generator smoke:
26|26|26|
27|27|27|```bash
28|28|28|./build/generate_bose build/UPC_CrossSection_AuAu_200GeV_XnXn_y1_dp001_nP200_Opt.root \
29|29|29|  review/logs/bose_weighted_1000_new.root 1000 123 1000 1
30|30|30|```
31|31|31|
32|32|32|Observed:
33|33|33|- Generated 1000/1000 events.
34|34|34|- Rejected candidates: 0.
35|35|35|- Wall time: about 0.95 s in this WSL run.
36|36|36|- 1000-event weight stats from an earlier baseline: w_min=4.1e-6, w_mean=0.982, w_rms=4.40, w_max=84.8, Neff~50/1000. The large weight tail is the key challenge for true unweighted generation.
37|37|37|
38|38|38|## b-integration convergence smoke
39|39|39|Command family:
40|40|40|
41|41|41|```bash
42|42|42|./build/generate_grid_opt AuAu 200 XnXn review/logs/grid_AuAu200_bstep${step}.root \
43|43|43|  6 16 20 ${step} 0.5 4 40
44|44|44|```
45|45|45|
46|46|46|Results:
47|47|47|
48|48|48|| bTailStep fm | sigma_y_window mb | momentum integral mb | wall s |
49|49|49||---:|---:|---:|---:|
50|50|50|| 8 | 0.764709 | 0.764709 | 0.78 |
51|51|51|| 4 | 0.764728 | 0.764728 | 0.88 |
52|52|52|| 2 | 0.764731 | 0.764731 | 1.03 |
53|53|53|
54|54|54|Conclusion for this AuAu200 XnXn smoke grid: tail step 4 fm is numerically stable at the few 1e-5 relative level versus 2 fm. This is not a full production validation, but it directly addresses the b-tail step concern for this configuration.
55|55|55|
56|56|56|## Other-energy smoke tests
57|57|57|Small grids ran for:
58|58|58|- AuAu 200 XnXn: loaded `UPC_Probabilities_AuAu_200GeV_Glauber_Final.root`.
59|59|59|- PbPb 5360 NoTag: loaded `UPC_Probabilities_PbPb_5360GeV_Glauber_Final.root`.
60|60|60|- PbPb 5020 NoTag: did not find a matching probability ROOT and fell back to code-level survival approximation.
61|61|61|
62|62|62|Important risk: arbitrary energies may silently use approximate fallback probabilities. For tagged triggers this can become invalid or throw if EMD histograms are missing. Production scans should either provide matching `UPC_Probabilities_*` files or explicitly document fallback usage.
63|63|63|
64|64|64|## Code changes made
65|65|65|Implemented a robust unweighted-event mode in `generate_bose`:
66|66|66|
67|67|67|```bash
68|68|68|generate_bose input.root output.root nEvents seed \
69|69|69|  [decay_norm_trials_per_mass] [bose_symmetrize=1] \
70|70|70|  [unweighted=0] [unweighting_trials=200000] [safety=1.25]
71|71|71|```
72|72|72|
73|73|73|Behavior:
74|74|74|- Default remains old weighted mode; existing commands are backward-compatible.
75|75|75|- If `unweighted=1`, the generator estimates a global accept-reject envelope from trial candidates.
76|76|76|- Accepted events are written with `event_weight=1`.
77|77|77|- ROOT metadata now records `unweighted_event_mode`, `unweighting_envelope`, `unweighting_trials`, `unweighting_safety_factor`, and `overweight_candidates`.
78|78|78|- If a later candidate exceeds the estimated envelope, generation aborts instead of silently clipping. This protects physics reliability; clipping would bias the unweighted sample.
79|79|79|
80|80|80|Unweighted smoke:
81|81|81|
82|82|82|```bash
83|83|83|./build/generate_bose build/UPC_CrossSection_AuAu_200GeV_XnXn_y1_dp001_nP200_Opt.root \
84|84|84|  review/logs/bose_unweighted_50_global50.root 50 123 1000 1 1 50000 50.0
85|85|85|```
86|86|86|
87|87|87|Observed:
88|88|88|- Generated 50/50 events.
89|89|89|- Rejected candidates: 584333.
90|90|90|- Overweight candidates: 0.
91|91|91|- All stored `event_weight` values were 1.
92|92|92|- Wall time: about 1.94 s.
93|93|93|
94|94|94|Interpretation: true unweighting is currently reliable but very inefficient because the decay/spin weight has a long tail. Weighted mode remains much more efficient statistically for high-statistics analysis unless a better proposal distribution is introduced.
95|95|95|
96|96|96|## Experimental / previous-result benchmark pass
97|97|97|
98|98|98|### STAR 4-prong observation
99|99|99|A literature check found STAR proceedings arXiv:0712.3527, which says the 2004 200 GeV Au-Au UPC data contain O(100) coherent 4-track candidates and a broad invariant-mass peak around 1500 MeV, possibly consistent with rho*(1450). The proceedings do not quote a final absolute 4-pion cross section, so it is a qualitative benchmark only.
100|100|100|
101|101|101|### STAR rho0 normalization cross-check
102|102|102|The STAR coherent rho0 AuAu 200 GeV XnXn measurement quoted in nucl-ex/0210028 is:
103|103|103|
104|104|104|- sigma(AuAu -> Au*_Xn Au*_Xn rho0) = 39.7 +/- 2.8 +/- 9.7 mb.
105|105|105|
106|106|106|Current/old rho-prime 4pi model outputs:
107|107|107|
108|108|108|- Current optimized grid in `build/UPC_CrossSection_AuAu_200GeV_XnXn_y1_dp001_nP200_Opt.root`: |y|<1 integral = 1.471 mb, d sigma/dy at y~0 = 0.759 mb.
109|109|109|- Older 2D model file `/mnt/d/Work/model_work/4pion/new-model/UPC_CrossSection_AuAu_200GeV_XnXn_2D.root`: full y integral = 2.890 mb, d sigma/dy at y~0 = 0.754 mb.
110|110|110|
111|111|111|Rough conclusion: predicted rho-prime->4pi XnXn is about 3.7% of STAR rho0 XnXn within |y|<1, or about 7.3% if the older full-y integral is used. That is not absurd for an excited rho/4pi channel and is qualitatively compatible with STAR only seeing O(100) 4-prong candidates. But because STAR does not quote a final 4pi cross section in this proceedings note, this is a sanity check, not a precision validation.
112|112|112|
113|113|113|## UPC probability ROOT/code recovery
114|114|114|Searched parent/sibling folders and found probability ROOT/code in:
115|115|115|
116|116|116|- `/mnt/d/Work/CodexCoding/doublerho_check/doublerho/UPC_Probabilities_PbPb_5020GeV_Glauber_Final.root`
117|117|117|- `/mnt/d/Work/CodexCoding/doublerho_check/doublerho/calc_UPC_EMD_PbPb_Glauber_new.C`
118|118|118|- `/mnt/d/Work/CodexCoding/doublerho_check/doublerho/calc_UPC_EMD_AuAu_Glauber_new.C`
119|119|119|
120|120|120|Actions:
121|121|121|
122|122|122|- Copied PbPb 5020 probability ROOT into this repo root: `UPC_Probabilities_PbPb_5020GeV_Glauber_Final.root`.
123|123|123|- Copied the two source macros into `review/probability_macros/` for provenance.
124|124|124|- Re-ran PbPb 5020 NoTag smoke. Before copying, the code fell back to approximate survival and gave sigma_y_window ~= 37.70 mb for the small smoke grid. After copying the proper ROOT, it loaded the probability file and gave sigma_y_window ~= 33.43 mb. This is a 12-13% effect even for NoTag, so missing probability ROOTs are not a harmless detail.
125|125|125|
126|126|126|## Iteration 2: faster unit-weight generation
127|127|127|Added `unweighting_mode=2`, a fast weighted-resampling unit-weight mode:
128|128|128|
129|129|129|```bash
130|130|130|generate_bose input.root output.root nEvents seed \
131|131|131|  [decay_norm_trials_per_mass] [bose_symmetrize=1] \
132|132|132|  [unweighted=0] [unweighting_trials=200000] [safety=1.25] [unweighting_mode]
133|133|133|```
134|134|134|
135|135|135|Modes:
136|136|136|
137|137|137|- `unweighting_mode=1`: exact accept-reject. Reliable but very slow for this long-tailed weight distribution.
138|138|138|- `unweighting_mode=2`: generate a weighted candidate pool, then resample events from the pool with probability proportional to candidate weight and write `event_weight=1`. This is much faster but is a pooled/resampled unit-weight approximation; events are not statistically independent beyond the pool size.
139|139|139|
140|140|140|Smoke comparison:
141|141|141|
142|142|142|- Exact accept-reject: 50 unit-weight events with 50k envelope trials and safety 50 took about 1.94 s and rejected 584333 candidates.
143|143|143|- Resampled mode: 1000 unit-weight events from a 5000-candidate pool took about 1.74 s.
144|144|144|- Unit-weight validation passed: 1000/1000 entries had `event_weight=1`.
145|145|145|- Weighted 5000-event reference vs resampled 1000-event means were close at smoke level:
146|146|146|  - weighted: <M>=1.5147, <y>=-0.0296, <pt>=0.05616
147|147|147|  - resampled: <M>=1.5041, <y>=-0.0391, <pt>=0.05501
148|148|148|
149|149|149|## Main reliability risks found
150|150|150|1. Exact true unweighted generation has a severe long-tail envelope problem. Low safety factors abort or would bias if clipped. Current implementation chooses abort over bias.
151|151|151|2. Fast resampling is much faster and practical, but it is a pooled approximation; production should use a large pool and compare key histograms to weighted reference.
152|152|152|3. Other collision energies can fall back to approximate survival treatment if the matching UPC probability ROOT is absent. PbPb 5020 smoke changed by ~12-13% after copying the proper probability ROOT.
153|153|153|4. `generate_grid_opt` uses fixed local transverse grid constants (`kNLoc=80`, `kLLocFm=40`) not exposed in the CLI. More convergence scans in local box/grid are needed for final production.
154|154|154|5. The optimized grid uses a 50-point omega LUT. Need a direct `generate_grid` vs `generate_grid_opt` comparison at representative settings to quantify interpolation bias.
155|155|155|6. Event analysis of Bose events relies on best rho/sigma pair assignment, not truth labels, so reconstructed angular modulations need systematic comparison to non-Bose/truth-labeled control.
156|156|156|
157|157|157|## Next recommended work
158|158|158|1. Add automated histogram comparison between weighted reference and resampled unit-weight output: M, y, pt, phi, cos2phi, pair masses.
159|159|159|2. Make probability-file fallback explicit in output metadata and optionally fatal for production mode.
160|160|160|3. Expose local grid size/box or at least record them as ROOT metadata.
161|161|161|4. Add automated convergence scans: bTailStep, bTailMax, nP/global_box, omega LUT points, and local grid size.
162|162|162|5. Improve resampling with stratified pools in M/y/pt cells to reduce duplicate events and long-tail variance.
163|163|163|
164|164|
165|165|## Iteration 3: validation tooling and production-safety metadata
166|166|
167|167|### Event-sample comparison tool
168|168|Added ROOT macro:
169|169|
170|170|```bash
171|171|root -l -b -q 'review/scripts/compare_event_samples.C("weighted.root","unit.root")'
172|172|```
173|173|
174|174|It compares a weighted reference TTree to a unit-weight/resampled TTree for:
175|175|
176|176|- `M_rhoprime`
177|177|- `y_rhoprime`
178|178|- `pt_rhoprime`
179|179|
180|180|It prints weighted/reference means, unit-sample means, RMS values, and shape chi2/ndf after histogram normalization.
181|181|
182|182|Smoke run:
183|183|
184|184|```bash
185|185|./build/generate_bose ... review/logs/compare_weighted_10000.root 10000 2026 1000 1
186|186|./build/generate_bose ... review/logs/compare_resample_5000.root 5000 2026 1000 1 1 20000 1.0 2
187|187|root -l -b -q 'review/scripts/compare_event_samples.C("review/logs/compare_weighted_10000.root","review/logs/compare_resample_5000.root")'
188|188|```
189|189|
190|190|Output summary:
191|191|
192|192|| observable | weighted mean | resampled mean | shape chi2/ndf |
193|193||---|---:|---:|---:|
194|194|| M_rhoprime | 1.53378 | 1.52560 | 1.23 |
195|195|| y_rhoprime | 0.02758 | 0.04248 | 0.68 |
196|196|| pt_rhoprime | 0.05329 | 0.05511 | 1.36 |
197|197|
198|198|Conclusion: for this smoke test, resampling reproduces the weighted shape at roughly chi2/ndf ~ O(1). This is a useful automated guard, not yet a final high-statistics validation.
199|199|
200|200|### Probability-file metadata and fatal production switch
201|201|`generate_grid_opt` now writes ROOT metadata:
202|202|
203|203|- `probability_file_loaded`
204|204|- `probability_fallback_used`
205|205|- `require_probability_file`
206|206|- `probability_file_path`
207|207|- `local_grid_n`
208|208|- `local_grid_box_fm`
209|209|- `local_grid_dx_fm`
210|210|- `omega_lut_points`
211|211|
212|212|The CLI now accepts a final optional argument:
213|213|
214|214|```bash
215|215|generate_grid_opt AuAu|PbPb sqrt_s_NN_GeV tag [output] [nM] [nP] [bCoreSteps] \
216|216|  [bTailStep] [yHalfWidth] [nY] [global_box_fm] [require_probability_file=0]
217|217|```
218|218|
219|219|Validation:
220|220|
221|221|```bash
222|222|./build/generate_grid_opt PbPb 5020 NoTag review/logs/grid_PbPb5020_metadata.root \
223|223|  4 10 10 8 0.5 4 40 1
224|224|```
225|225|
226|226|Metadata readback:
227|227|
228|228|```text
229|229|probability_file_loaded=1
230|230|probability_fallback_used=0
231|231|require_probability_file=1
232|232|local_grid_n=80
233|233|omega_lut_points=50
234|234|probability_file_path=UPC_Probabilities_PbPb_5020GeV_Glauber_Final.root
235|235|```
236|236|
237|237|Fatal-missing-probability check:
238|238|
239|239|```bash
240|240|./build/generate_grid_opt PbPb 1234 NoTag review/logs/grid_missingprob_shouldfail.root \
241|241|  2 4 2 20 0.1 2 40 1
242|242|```
243|243|
244|244|Output:
245|245|
246|246|```text
247|247|generate_grid_opt: cannot load required UPC probability ROOT for PbPb @ 1234.000000
248|248|exit_code=1
249|249|```
250|250|
251|251|So production jobs can now fail fast instead of silently using fallback probabilities.
252|252|
253|253|### Direct vs optimized grid check
254|254|Small direct-vs-optimized comparison:
255|255|
256|256|```bash
257|257|./build/generate_grid     AuAu 200 XnXn review/logs/generate_grid_AuAu_cmp.root     3 6 4 20 0.2 2 40
258|258|./build/generate_grid_opt AuAu 200 XnXn review/logs/generate_grid_opt_AuAu_cmp.root 3 6 4 20 0.2 2 40
259|259|```
260|260|
261|261|Results:
262|262|
263|263|| implementation | sigma_y_window mb | dsigma/dy y~0 mb | spatial mb | momentum mb | wall s |
264|264||---|---:|---:|---:|---:|---:|
265|265|| direct | 0.369790 | 0.924474 | 0.376577 | 0.369790 | 0.35 |
266|266|| opt | 0.369757 | 0.924392 | 0.376544 | 0.369757 | 0.32 |
267|267|
268|268|Relative difference in integrated momentum result is about 9e-5 for this smoke configuration. The omega-LUT optimized path is therefore numerically consistent with direct evaluation at this small representative point.
269|269|
270|
271|## Iteration 4: event-generator efficiency and duplicate control
272|
273|### Baseline profiling
274|Commands were timed with `/usr/bin/time` using the AuAu200 XnXn optimized input grid.
275|
276|| mode | command meaning | output events | pool/trials | wall s | max RSS |
277||---|---|---:|---:|---:|---:|
278|| weighted | ordinary weighted generation | 20000 | n/a | 7.36 | 556 MB |
279|| mode 2 | resample with replacement | 20000 | 50000 | 1.22 | 582 MB |
280|| mode 3 | weighted reservoir without replacement | 20000 | 50000 | 1.18 | 578 MB |
281|
282|The weighted TTree output path is relatively slow for large unique-event samples. Unit-weight pool/reservoir modes are much faster because the expensive candidate generation is capped by the pool size and, in mode 2, the output events can be sampled from the pool.
283|
284|### New unit-weight mode 3: weighted reservoir
285|Added:
286|
287|```bash
288|unweighting_mode=3
289|```
290|
291|This performs one pass over `unweighting_trials` weighted candidates and keeps a weighted reservoir of size `nEvents` using the Efraimidis-Spirakis key `log(U)/w`. It writes unit-weight events without replacement from the selected reservoir.
292|
293|Mode comparison:
294|
295|- `mode=2`: fastest practical resampling with replacement; can duplicate events if the weight distribution has a long tail.
296|- `mode=3`: similar speed, bounded memory O(nEvents), and no duplicate selected candidates in the tested kinematic key. Better for analyses sensitive to repeated identical events.
297|
298|Duplicate check for 20k output events from a 50k pool/trial budget:
299|
300|| mode | unique kinematic keys | duplicate fraction |
301||---|---:|---:|
302|| mode 2 resampling | 9317 / 20000 | 0.534 |
303|| mode 3 reservoir | 20000 / 20000 | 0.000 |
304|
305|Here the key was `(M_rhoprime, y_rhoprime, pt_rhoprime, m_rho, m_sigma)` rounded to printed double precision. This is not a full event identity hash, but it is enough to show that mode 2 heavily repeats high-weight candidates while mode 3 avoids that pathology.
306|
307|### Distribution self-check against weighted reference
308|Weighted reference: 20k ordinary weighted events. Unit-weight tests: 20k events with 50k pool/trials.
309|
310|| observable | mode 2 chi2/ndf | mode 3 chi2/ndf | comment |
311||---|---:|---:|---|
312|| M_rhoprime | 1.56 | 2.83 | both usable for smoke; mode 3 mass shape slightly worse at this pool size |
313|| y_rhoprime | 1.71 | 0.87 | both OK |
314|| pt_rhoprime | 2.34 | 1.73 | both OK |
315|
316|Interpretation: mode 2 best approximates the exact weighted distribution for a fixed small pool because it samples with replacement proportional to weights, but repeats many events. Mode 3 trades a little shape precision for zero duplicates and similar speed. Production recommendation: use mode 2 for fast shape prototyping; use mode 3 or larger-pool mode 2 for analyses where duplicate events matter.
317|
318|### Next generator-optimization target
319|The remaining bottleneck is not just micro-optimization; it is proposal quality. The decay/spin weight has a long tail. The next meaningful speed/reliability improvement is stratified pooling by M/y/pt or adaptive pool enlargement until weighted-vs-unit chi2 is acceptable in all monitored observables.
320|

## Iteration 5: ROOT output compression and TTree write speed

Added CLI/config parameter:

```bash
[compression_level=1]
```

to `generate_bose`, recorded as ROOT metadata `output_compression_level`. The generator now explicitly calls `TFile::SetCompressionLevel(cfg.output_compression_level)` after opening the output file.

### Speed / size test
AuAu200 XnXn input, 20k output events:

| mode | compression | wall s | output size |
|---|---:|---:|---:|
| weighted | 0 | 1.07 | 10 MB |
| weighted | 1 | 1.10 | 5.5 MB |
| weighted | 5 | 1.09 | 5.5 MB |
| reservoir mode 3 | 0 | 1.27 | 10 MB |
| reservoir mode 3 | 1 | 1.20 | 5.4 MB |

A previous weighted 20k profile before explicit compression control showed a much slower wall time (~7.36 s). After explicitly setting the compression level, weighted 20k is ~1.1 s. For these TLorentzVector-heavy TTrees, compression level 1 keeps file size about half of uncompressed with no measurable speed penalty in this WSL smoke test.

### Current event-generator recommendation
- For exact weighted studies: use default weighted mode with compression level 1.
- For fast unit-weight prototyping: `unweighted=1, unweighting_mode=2` with a sufficiently large pool.
- For unit-weight samples where duplicate events are problematic: `unweighted=1, unweighting_mode=3`.
- Use `review/scripts/compare_event_samples.C` against a weighted reference after changing pool size or mode.

Example fast no-duplicate unit-weight command:

```bash
./build/generate_bose input.root output.root 20000 6161 1000 1 1 50000 1.0 3 1
```
