#manifest
echo -e "sample-id\tabsolute-filepath" > manifest
for i in *fastq.gz; do 
	path=$(readlink -f $i); id=$(basename $i | cut -d '.' -f 1)
    echo -e $id'\t'$path >>manifest
done

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# demux-view
qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization demux.qzv

#deblur
qiime quality-filter q-score \
  --i-demux single-end-demux.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza

qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv

qiime demux summarize \
  --i-data demux-filtered.qza \
  --o-visualization demux-filtered.qzv

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 180 \
  --p-jobs-to-start 0 \
  --p-no-hashed-feature-ids \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza

qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

mv rep-seqs-deblur.qza rep-seqs.qza

mv table-deblur.qza table.qza

#FeatureTable and FeatureData summaries
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv

qiime feature-table tabulate-seqs \
 --i-data rep-seqs.qza \
 --o-visualization rep-seqs.qzv

#Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --p-n-threads 28 \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 37000 \
  --m-metadata-file metadata.txt \
  --output-dir core-metrics-results

qiime feature-classifier classify-sklearn \
  --i-classifier gg2-v3v4-classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 28 \
  --o-classification taxonomy-new.qza
