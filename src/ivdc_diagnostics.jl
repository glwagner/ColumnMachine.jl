using ClimateMachine.Writers

function ivdc_default_diagnostics_init()
end
function ivdc_default_diagnostics_fini()
end
function ivdc_default_diagnostics_collect()
    mpicomm = Diagnostics.Settings.mpicomm
    dg1 = Diagnostics.Settings.dg.modeldata.diagsconf.myivdcdg
    dg2 = Diagnostics.Settings.dg.modeldata.diagsconf.myivdcdg2
    Q = Diagnostics.Settings.dg.modeldata.diagsconf.Q
end

function setup_single_column_default_diagnostics(::IVDCConfigType, interval, out_pref; writer=NetCDFWriter(),interpol=nothing)

	return DiagnosticsGroup(
				"IVDCdefault",
				ivdc_default_diagnostics_init,
				ivdc_default_diagnostics_fini,
				ivdc_default_diagnostics_collect,
				interval,
				out_pref,
				writer,
				interpol,
				)

end
