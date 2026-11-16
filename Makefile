all := libareole
LDLIBS-libareole := -lm

WASM_SRCS := src/areole.c npm/wasm/areole_wasm.c

include ../mk/include.mk
LDLIBS-alpine :=

CFLAGS := -I./include -I../libislet/include -O2 -Wno-overlength-strings ${CFLAGS} $(EXTRA_CFLAGS)

sync-npm-headers: include/ttypt/areole.h include/ttypt/areole_glsl.h
	@mkdir -p npm/include/ttypt
	@cp -f include/ttypt/areole.h include/ttypt/areole_glsl.h npm/include/ttypt/

all: sync-npm-headers

wasm: sync-npm-headers ${wasmtarget}
	@cp -f npm/dist/areole.wasm npm/dist/chunk.wasm
	@if [ -d ../libchunk/npm/dist ]; then cp -f npm/dist/areole.wasm ../libchunk/npm/dist/chunk.wasm; fi

test: all bin/test_areole${EXE} bin/test_glsl${EXE}
	./bin/test_areole${EXE}
	./bin/test_glsl${EXE}

bin/test_areole${EXE}: tests/test_areole.c src/areole.c include/ttypt/areole.h
	@mkdir -p bin
	${cc} ${CFLAGS} -o $@ tests/test_areole.c src/areole.c -lm

bin/test_glsl${EXE}: tests/test_glsl.c include/ttypt/areole_glsl.h
	@mkdir -p bin
	${cc} ${CFLAGS} -o $@ tests/test_glsl.c -lm

${DESTDIR}${PREFIX}/lib/pkgconfig/areole.pc: areole.pc
	install -d ${DESTDIR}${PREFIX}/lib/pkgconfig
	install -m 644 areole.pc $@

install: ${DESTDIR}${PREFIX}/lib/pkgconfig/areole.pc

NPM ?= npm

npm-build: wasm
	@if command -v $(NPM) >/dev/null 2>&1; then \
		(cd npm && $(NPM) run build 2>/dev/null || true); \
	fi

npm-test: npm-build
	@if command -v $(NPM) >/dev/null 2>&1; then \
		(cd npm && $(NPM) test); \
	fi

clean-extra::
	rm -rf bin npm/build npm/dist/*.wasm

objects-set.mk: Makefile

.PHONY: all sync-npm-headers test wasm npm-build npm-test install
