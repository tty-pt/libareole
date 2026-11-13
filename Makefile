# libareole Makefile — Generic Distance-Adaptive Multi-Scale Spatial Quadtree & LOD Engine
PREFIX ?= /usr/local
DESTDIR ?=

CC ?= gcc
CFLAGS ?= -Wall -Wextra -Werror -O2 -Iinclude -I../libislet/include
LDLIBS ?= -lm

WASM_CC ?= clang
WASM_CFLAGS ?= --target=wasm32-wasi -ffreestanding -flto -O3 -D__wasm__ -fvisibility=hidden -I./include -I../libislet/include -I../libstele/include -I../mk/include -I/usr/include
WASM_LDFLAGS ?= -mexec-model=reactor -Wl,--gc-sections -Wl,--no-entry -Wl,--allow-undefined \
	-Wl,--import-memory -Wl,-z,stack-size=32768 -Wl,--global-base=1769472 \
	-Wl,--wrap=malloc -Wl,--wrap=free -Wl,--wrap=calloc -Wl,--wrap=realloc \
	-Wl,--export=ar_wasm_create \
	-Wl,--export=ar_wasm_destroy \
	-Wl,--export=ar_wasm_clear \
	-Wl,--export=ar_wasm_update \
	-Wl,--export=ar_wasm_probe \
	-Wl,--export=ar_wasm_get_node_ptr \
	-Wl,--export=ar_wasm_get_node_count

NPM ?= npm

all: lib/libareole.so test wasm

lib/libareole.so: src/areole.c include/ttypt/areole.h
	@mkdir -p lib
	$(CC) $(CFLAGS) -fPIC -shared -o $@ src/areole.c $(LDLIBS)

test: bin/test_areole bin/test_glsl
	./bin/test_areole
	./bin/test_glsl

bin/test_areole: tests/test_areole.c src/areole.c include/ttypt/areole.h
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ tests/test_areole.c src/areole.c $(LDLIBS)

bin/test_glsl: tests/test_glsl.c include/ttypt/areole_glsl.h
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ tests/test_glsl.c $(LDLIBS)

wasm: npm/dist/areole.wasm

npm/dist/areole.wasm: npm/wasm/areole_wasm.c src/areole.c include/ttypt/areole.h ../mk/wasm_alloc_shim.c
	@mkdir -p npm/dist
	$(WASM_CC) $(WASM_CFLAGS) $(WASM_LDFLAGS) npm/wasm/areole_wasm.c src/areole.c ../mk/wasm_alloc_shim.c -o $@
	@ln -sf areole.wasm npm/dist/chunk.wasm
	@if [ -d ../libchunk/npm/dist ]; then cp -f $@ ../libchunk/npm/dist/chunk.wasm; fi

npm-build: wasm
	@if command -v $(NPM) >/dev/null 2>&1; then \
		(cd npm && $(NPM) run build 2>/dev/null || true); \
	fi

npm-test: npm-build
	@if command -v $(NPM) >/dev/null 2>&1; then \
		(cd npm && node --test tests/areole.test.mjs); \
	fi

install: lib/libareole.so areole.pc
	install -d $(DESTDIR)$(PREFIX)/include/ttypt
	install -m 644 include/ttypt/areole.h $(DESTDIR)$(PREFIX)/include/ttypt/
	install -m 644 include/ttypt/areole_glsl.h $(DESTDIR)$(PREFIX)/include/ttypt/
	install -m 644 $(DESTDIR)$(PREFIX)/include/ttypt/
	install -m 644 $(DESTDIR)$(PREFIX)/include/ttypt/
	install -d $(DESTDIR)$(PREFIX)/lib
	install -m 755 lib/libareole.so $(DESTDIR)$(PREFIX)/lib/
	install -d $(DESTDIR)$(PREFIX)/lib/pkgconfig
	install -m 644 areole.pc $(DESTDIR)$(PREFIX)/lib/pkgconfig/

clean:
	rm -rf bin lib npm/dist/*.wasm

.PHONY: all test wasm npm-build npm-test install clean
