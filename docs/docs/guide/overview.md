# Overview

A vAIRRify audit follows three stages: **prepare the data**, **define the study design**, and run **unsupervised** and **supervised** confounder analyses on the locked sample set. Click any box in the diagram to jump to its explanation.

<div class="vairrify-workflow">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 1032" style="width: 100%; height: auto; display: block;" role="img" aria-label="vAIRRify step-by-step workflow">
  <defs>
    <style>
      .vf-bg { fill: #ffffff; }
      .vf-font { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
      .vf-mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
      .vf-group-prep { fill: #f6f8fa; stroke: #8c959f; stroke-width: 2.4; rx: 10; }
      .vf-group-design { fill: #f0f6fc; stroke: #2a78d6; stroke-width: 2.4; rx: 10; }
      .vf-group-unsup { fill: #f3f8f7; stroke: #3c9b85; stroke-width: 2.4; rx: 10; }
      .vf-group-sup { fill: #fdf8f2; stroke: #c9863b; stroke-width: 2.4; rx: 10; }
      .vf-group-label { font-size: 11.5px; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase; fill: #57606a; }
      .vf-group-label-design { fill: #164e96; }
      .vf-group-label-unsup { fill: #1c6b5a; }
      .vf-group-label-sup { fill: #9a5514; }
      .vf-node { fill: #ffffff; stroke: #8c959f; stroke-width: 1.3; rx: 8; }
      .vf-node-opt { fill: #ffffff; stroke: #8c959f; stroke-width: 1.3; stroke-dasharray: 5 3; rx: 8; }
      .vf-node-design { fill: #ffffff; stroke: #2a78d6; stroke-width: 1.35; rx: 8; }
      .vf-node-unsup { fill: #ffffff; stroke: #2d947c; stroke-width: 1.35; rx: 8; }
      .vf-node-sup { fill: #ffffff; stroke: #d47b27; stroke-width: 1.35; rx: 8; }
      .vf-node-footer { fill: #f6f8fa; stroke: #8c959f; stroke-width: 1.3; rx: 7; }
      .vf-role-pill { fill: #f5f9fe; stroke: #bdd6f2; stroke-width: 1; rx: 5; }
      .vf-role-pill-opt { fill: #f5f9fe; stroke: #bdd6f2; stroke-width: 1; stroke-dasharray: 4 2.5; rx: 5; }
      .vf-role-name { font-size: 11px; font-weight: 700; fill: #164e96; }
      .vf-role-desc { font-size: 10px; fill: #57606a; }
      .vairrify-workflow a { text-decoration: none; cursor: pointer; }
      .vairrify-workflow a:hover rect { filter: brightness(0.96); stroke-width: 2; }
      .vf-title { font-size: 14px; font-weight: 700; fill: #1f2328; }
      .vf-code { font-size: 11.5px; fill: #2a5ea8; }
      .vf-code-sm { font-size: 10.8px; fill: #2a5ea8; }
      .vf-sub { font-size: 11.5px; fill: #57606a; }
      .vf-badge-opt { fill: #eef1f4; stroke: #d0d7de; stroke-width: 0.8; rx: 7; }
      .vf-badge-opt-text { font-size: 10px; font-weight: 600; fill: #57606a; }
      .vf-edge-label { font-size: 11px; font-weight: 600; fill: #57606a; }
      .vf-edge-label-bg-sup { fill: #fdf8f2; rx: 3; }
      .vf-arrow { fill: none; stroke: #57606a; stroke-width: 1.6; stroke-linecap: round; stroke-linejoin: round; }
      .vf-arrow-sup { fill: none; stroke: #c46d1b; stroke-width: 1.6; stroke-linecap: round; stroke-linejoin: round; }
      .vf-stage-arrow { fill: #cbd2d9; }
    </style>
    <marker id="vf-ah-gray" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 1 1.5 L 9 5 L 1 8.5 z" fill="#57606a" />
    </marker>
    <marker id="vf-ah-orange" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 1 1.5 L 9 5 L 1 8.5 z" fill="#c46d1b" />
    </marker>
  </defs>
  <rect class="vf-bg" x="0" y="0" width="680" height="1032" />
  <rect class="vf-group-prep" x="25" y="16" width="630" height="240" />
  <text class="vf-font vf-group-label" x="41" y="37">1. Prepare data</text>
  <a href="../../getting-started/">
    <rect class="vf-node" x="41" y="50" width="186" height="78" />
    <text class="vf-font vf-title" x="134" y="72" text-anchor="middle">Ingest</text>
    <text class="vf-mono vf-code" x="134" y="90" text-anchor="middle">tcrio · DatasetIngester</text>
    <text class="vf-font vf-sub" x="134" y="106" text-anchor="middle">standardise &amp; filter</text>
    <text class="vf-font vf-sub" x="134" y="120" text-anchor="middle">clonotype tables</text>
  </a>
  <a href="../../getting-started/">
    <rect class="vf-node" x="247" y="50" width="186" height="78" />
    <text class="vf-font vf-title" x="340" y="74" text-anchor="middle">Add metadata</text>
    <text class="vf-mono vf-code" x="340" y="94" text-anchor="middle">patient &amp; repertoire</text>
    <text class="vf-font vf-sub" x="340" y="113" text-anchor="middle">overview() · with_metadata()</text>
  </a>
  <a href="../../getting-started/">
    <rect class="vf-node" x="453" y="50" width="186" height="78" />
    <text class="vf-font vf-title" x="546" y="74" text-anchor="middle">Select locus</text>
    <text class="vf-mono vf-code" x="546" y="94" text-anchor="middle">TRA · TRB · IGH · ...</text>
    <text class="vf-font vf-sub" x="546" y="113" text-anchor="middle">one locus per audit</text>
  </a>
  <line class="vf-arrow" x1="227" y1="89" x2="243" y2="89" marker-end="url(#vf-ah-gray)" />
  <line class="vf-arrow" x1="433" y1="89" x2="449" y2="89" marker-end="url(#vf-ah-gray)" />
  <path class="vf-arrow" d="M 546 128 L 546 144 L 134 144 L 134 156" marker-end="url(#vf-ah-gray)" />
  <a href="../../getting-started/">
    <rect class="vf-node-opt" x="41" y="160" width="186" height="80" />
    <rect class="vf-badge-opt" x="155" y="152" width="60" height="16" />
    <text class="vf-font vf-badge-opt-text" x="185" y="164" text-anchor="middle">optional</text>
    <text class="vf-font vf-title" x="134" y="186" text-anchor="middle">Check overlap</text>
    <text class="vf-mono vf-code" x="134" y="205" text-anchor="middle">OverlapAnalyzer</text>
    <text class="vf-font vf-sub" x="134" y="224" text-anchor="middle">sample contamination</text>
  </a>
  <a href="../../feature-blocks/">
    <rect class="vf-node" x="247" y="160" width="186" height="80" />
    <text class="vf-font vf-title" x="340" y="180" text-anchor="middle">Compute repertoire</text>
    <text class="vf-font vf-title" x="340" y="196" text-anchor="middle">properties</text>
    <text class="vf-mono vf-code" x="340" y="213" text-anchor="middle">reps = vf.load(...)</text>
    <text class="vf-font vf-sub" x="340" y="229" text-anchor="middle">V/J, VJ bias, VDJ, depth, HLA</text>
  </a>
  <a href="../qc/">
    <rect class="vf-node" x="453" y="160" width="186" height="80" />
    <text class="vf-font vf-title" x="546" y="182" text-anchor="middle">Quality control (QC)</text>
    <text class="vf-mono vf-code" x="546" y="200" text-anchor="middle">qc = reps.qc()</text>
    <text class="vf-font vf-sub" x="546" y="217" text-anchor="middle">low depth · missing</text>
    <text class="vf-font vf-sub" x="546" y="231" text-anchor="middle">genes · outliers</text>
  </a>
  <line class="vf-arrow" x1="227" y1="200" x2="243" y2="200" marker-end="url(#vf-ah-gray)" />
  <line class="vf-arrow" x1="433" y1="200" x2="449" y2="200" marker-end="url(#vf-ah-gray)" />
  <polygon class="vf-stage-arrow" points="333,260 347,260 347,274 358,274 340,289 322,274 333,274" />
  <text class="vf-mono vf-edge-label" x="366" y="278">exclude=qc.failed</text>
  <rect class="vf-group-design" x="25" y="293" width="630" height="176" />
  <text class="vf-font vf-group-label vf-group-label-design" x="41" y="314">2. Define study design</text>
  <a href="../design/">
    <rect class="vf-node-design" x="41" y="327" width="598" height="128" />
    <text class="vf-mono vf-code-sm" x="340" y="347" text-anchor="middle">design = reps.design(label=..., batch=..., other_technical=[...], covariates=[...], controls=[...])</text>
    <rect class="vf-role-pill" x="55" y="358" width="180" height="34" />
    <text class="vf-mono vf-role-name" x="145" y="372" text-anchor="middle">label</text>
    <text class="vf-font vf-role-desc" x="145" y="386" text-anchor="middle">target to predict (e.g. disease)</text>
    <rect class="vf-role-pill" x="250" y="358" width="180" height="34" />
    <text class="vf-mono vf-role-name" x="340" y="372" text-anchor="middle">batch</text>
    <text class="vf-font vf-role-desc" x="340" y="386" text-anchor="middle">main technical grouping</text>
    <rect class="vf-role-pill" x="445" y="358" width="180" height="34" />
    <text class="vf-mono vf-role-name" x="535" y="372" text-anchor="middle">other_technical</text>
    <text class="vf-font vf-role-desc" x="535" y="386" text-anchor="middle">depth, run, experiment</text>
    <rect class="vf-role-pill" x="55" y="400" width="180" height="34" />
    <text class="vf-mono vf-role-name" x="145" y="414" text-anchor="middle">covariates</text>
    <text class="vf-font vf-role-desc" x="145" y="428" text-anchor="middle">age, sex, ancestry</text>
    <rect class="vf-role-pill" x="250" y="400" width="180" height="34" />
    <text class="vf-mono vf-role-name" x="340" y="414" text-anchor="middle">controls</text>
    <text class="vf-font vf-role-desc" x="340" y="428" text-anchor="middle">negative controls (random_control)</text>
    <rect class="vf-role-pill-opt" x="445" y="400" width="180" height="34" />
    <text class="vf-mono vf-role-name" x="535" y="414" text-anchor="middle">cv_column</text>
    <text class="vf-font vf-role-desc" x="535" y="428" text-anchor="middle">cross-validation folds (optional)</text>
    <text class="vf-font vf-sub" x="340" y="448" text-anchor="middle">Excludes failed QC samples, merges rare covariate levels, and locks one sample set for Step 3</text>
  </a>
  <polygon class="vf-stage-arrow" points="333,473 347,473 347,487 358,487 340,502 322,487 333,487" />
  <rect class="vf-group-unsup" x="25" y="506" width="630" height="148" />
  <text class="vf-font vf-group-label vf-group-label-unsup" x="41" y="527">3a. Unsupervised confounder analysis</text>
  <a href="../balance/">
    <rect class="vf-node-unsup" x="41" y="540" width="190" height="98" />
    <text class="vf-font vf-title" x="136" y="564" text-anchor="middle">Balance over batches</text>
    <text class="vf-mono vf-code" x="136" y="584" text-anchor="middle">design.balance()</text>
    <text class="vf-font vf-sub" x="136" y="606" text-anchor="middle">Cramér's V · single- vs.</text>
    <text class="vf-font vf-sub" x="136" y="622" text-anchor="middle">mixed-label batches</text>
  </a>
  <a href="../pca/">
    <rect class="vf-node-unsup" x="245" y="540" width="190" height="98" />
    <text class="vf-font vf-title" x="340" y="564" text-anchor="middle">PCA &amp; loadings</text>
    <text class="vf-mono vf-code" x="340" y="584" text-anchor="middle">design.pca()</text>
    <text class="vf-font vf-sub" x="340" y="606" text-anchor="middle">Batch vs. label clusters ·</text>
    <text class="vf-font vf-sub" x="340" y="622" text-anchor="middle">block shares · germline V</text>
  </a>
  <a href="../variance/">
    <rect class="vf-node-unsup" x="449" y="540" width="190" height="98" />
    <text class="vf-font vf-title" x="544" y="564" text-anchor="middle">Sources of variation</text>
    <text class="vf-mono vf-code" x="544" y="584" text-anchor="middle">design.variance()</text>
    <text class="vf-font vf-sub" x="544" y="606" text-anchor="middle">PVCA share · PERMANOVA</text>
    <text class="vf-font vf-sub" x="544" y="622" text-anchor="middle">R² (alone vs. adjusted)</text>
  </a>
  <polygon class="vf-stage-arrow" points="333,658 347,658 347,672 358,672 340,687 322,672 333,672" />
  <rect class="vf-group-sup" x="25" y="691" width="630" height="270" />
  <text class="vf-font vf-group-label vf-group-label-sup" x="41" y="712">3b. Supervised confounder analysis</text>
  <a href="../baseline/">
    <rect class="vf-node-sup" x="41" y="725" width="598" height="86" />
    <text class="vf-font vf-title" x="340" y="749" text-anchor="middle">Classification based on general repertoire properties</text>
    <text class="vf-mono vf-code" x="340" y="768" text-anchor="middle">base = design.baseline()</text>
    <text class="vf-font vf-sub" x="340" y="786" text-anchor="middle">General · Predicted HLA · General + HLA</text>
    <text class="vf-font vf-sub" x="340" y="801" text-anchor="middle">Out-of-fold random-forest AUC &amp; AP (label, batch, controls)</text>
  </a>
  <path class="vf-arrow-sup" d="M 340 811 L 340 833 L 185 833 L 185 851" marker-end="url(#vf-ah-orange)" />
  <path class="vf-arrow-sup" d="M 340 833 L 495 833 L 495 851" marker-end="url(#vf-ah-orange)" />
  <rect class="vf-edge-label-bg-sup" x="346" y="814" width="34" height="16" />
  <text class="vf-mono vf-edge-label" x="363" y="826" text-anchor="middle">base</text>
  <a href="../within-batch/">
    <rect class="vf-node-sup" x="41" y="855" width="288" height="90" />
    <text class="vf-font vf-title" x="185" y="879" text-anchor="middle">Within-batch evaluation</text>
    <text class="vf-mono vf-code" x="185" y="898" text-anchor="middle">design.within_batch(base)</text>
    <text class="vf-font vf-sub" x="185" y="918" text-anchor="middle">Does baseline AUC hold inside</text>
    <text class="vf-font vf-sub" x="185" y="933" text-anchor="middle">mixed-label batches?</text>
  </a>
  <a href="../within-batch/">
    <rect class="vf-node-sup" x="351" y="855" width="288" height="90" />
    <text class="vf-font vf-title" x="495" y="879" text-anchor="middle">Feature attribution (SHAP)</text>
    <text class="vf-mono vf-code" x="495" y="898" text-anchor="middle">design.shap(base)</text>
    <text class="vf-font vf-sub" x="495" y="918" text-anchor="middle">Top features · value &amp;</text>
    <text class="vf-font vf-sub" x="495" y="933" text-anchor="middle">batch-coloured beeswarms</text>
  </a>
  <rect class="vf-node-footer" x="25" y="975" width="630" height="42" />
  <text class="vf-font vf-sub" x="340" y="1000" text-anchor="middle">Every step returns a table and plot (<tspan class="vf-mono vf-code">vf.plots.*</tspan>) and exports via <tspan class="vf-mono vf-code">vf.save(...)</tspan>. Click any box for its guide.</text>
</svg>
</div>
