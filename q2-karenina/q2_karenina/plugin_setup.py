import qiime.plugin

import q2_karenina
from q2_karenina import fit_timeseries
from q2_types.ordination import PCoAResults as pcoa
from qiime2.plugin import Metadata, Str, Choices


plugin = qiime.plugin.Plugin(
    name='karenina',
    version=q2_karenina.__version__,
    website='https://github.com/zaneveld/karenina',
    package='q2_karenina',
    user_support_text=None,
	description="This script simulates microbiome " +
    "change over time using Ornstein-Uhlenbeck (OU) models.  These are " +
    "similar to Brownian motion models, with the exception that they " +
    "include reversion to a mean. Output is a tab-delimited data table " +
    "and figures.",
    citation_text=None
)

plugin.visualizers.register_function(
    function=q2_karenina.fit_timeseries,
	inputs={
		'pcoa' : pcoa
    },
    parameters={
        'method' : Str % Choices({'basinhopping'})
		'metadata' : Metadata
		'individual' : Str
		'timepoint' : Str
		'treatment' : Str
    },
	parameter_descriptions = {
		'method' : 'global optimization method'
		'metadata' : 'Sample metadata'
		'individual' : 'individual column identifier'
		'timepoint' : 'timepoint column identifier'
		'treatment' : 'treatment column identifier'
	}
    name='Fit OU Models to PCoA Ordination output',
    description='This visualizer generates OU model parameters for PCoA output'
                'data, for each individual and each defined treatment cohort.'
)
