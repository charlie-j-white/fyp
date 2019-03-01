rm MERGED.f
#cat *.f > MERGED.f
find . -maxdepth 1 -iname '*.f' -not -name 'wrapper.f' -exec cat {} +>MERGED.f
