let NERDTreeIgnore=['\.pdf']


let g:makeprg="ninja"
"let g:cmake_opts="-G Ninja"

let s:build_dir=expand('<sfile>:p:h')."/build"
let g:dispatch="cd ".s:build_dir."; ninja test"

map <F9> :Dispatch<CR>

:CMake

call unite#custom#source('file_rec', 'ignore_pattern', '/build/.*' )
