project('opencalphad', ['fortran', 'c'], default_options: ['fortran_std=legacy'], version: '6.057', license: 'GPL-3.0', meson_version: '>=0.49')

lapack_dep = dependency('lapack')

buildtype = get_option('buildtype')
if buildtype == 'release'
  add_project_arguments(['-Ofast', '-march=native', '-pthread'], language : 'fortran')
endif

main_sources = files()
ocasi_sources = files()
tinyfiledialogs_sources = files()

subdir('src')
subdir('src/utilities')
subdir('src/utilities/GETKEY')
subdir('src/utilities/tinyfiledialogs')
subdir('src/minimizer')
subdir('src/models')
subdir('src/numlib')
subdir('src/stepmapplot')
subdir('src/userif')
subdir('OCisoCbinding')

ocasi_lib = shared_library('OCASI', ocasi_sources, dependencies: [lapack_dep])
ocasi_dep = declare_dependency(link_with: ocasi_lib)
tinyfiledialogs_lib = shared_library('tinyfiledialogs', tinyfiledialogs_sources)

executable('OC6', main_sources, dependencies : [ocasi_dep], link_with: [tinyfiledialogs_lib])

