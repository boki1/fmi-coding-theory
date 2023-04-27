# Solution to Assignment 1, Problem 2

*Implementation of linear codes*

-----------

**Building**

_Requires:_
- `conan, 1.57`
- `cmake, 3.22`
- `g++, 11.3.0`
- `ccmake`, (optional)

_Note:_
The versions of the tools mentioned are the ones that I have used during development. The implementation may work with others also, but that is not guaranteed.

```bash
$ output_dir=build
$ conan install . --output-folder=${output_dir} -if=${output_dir} --build=missing
$ cmake -S. -B${output_dir} -GNinja
$ ninja -C ${output_dir}

# If you want to build the demo also, you can change the `LCODES_BUILD_DEMO` option using `ccmake`.
$ ccmake ${output_dir}
```
-----------

**LICENSE**

MIT License, Copyright (c) 2023 Kristiyan Stoimenov
