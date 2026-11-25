all := libareole
LDLIBS-libareole := -lm

WASM_SRCS := src/areole.c npm/wasm/areole_wasm.c

include ../mk/include.mk
LDLIBS-alpine :=

CFLAGS := -I./include -I../libislet/include -O2 -Wno-overlength-strings ${CFLAGS} $(EXTRA_CFLAGS)

sync-npm-headers: include/ttypt/areole.h include/ttypt/areole_glsl.h
	@install -d npm/include/ttypt
	@install -m 644 include/ttypt/areole.h include/ttypt/areole_glsl.h npm/include/ttypt/

all: sync-npm-headers

wasm: sync-npm-headers ${wasmtarget}
	@install -m 644 npm/dist/areole.wasm npm/dist/chunk.wasm
	@if [ -d ../libchunk/npm/dist ]; then install -m 644 npm/dist/areole.wasm ../libchunk/npm/dist/chunk.wasm; fi

AREOLE_PART_SRCS := src/areole.c src/areole_internal.h src/ar_math.c src/ar_occ.c src/ar_pool.c src/ar_populate.c src/ar_solver.c src/ar_probe.c

test: all bin/test_areole${EXE} bin/test_glsl${EXE} bin/bench_areole${EXE}
	./bin/test_areole${EXE}
	./bin/test_glsl${EXE}
	./bin/bench_areole${EXE}

bench: all bin/bench_areole${EXE}
	./bin/bench_areole${EXE}

bin/test_areole${EXE}: tests/test_areole.c $(AREOLE_PART_SRCS) include/ttypt/areole.h
	@mkdir -p bin
	${cc} ${CFLAGS} -o $@ tests/test_areole.c src/areole.c ../libstele/src/stele_wq.c -lm -lpthread

bin/bench_areole${EXE}: tests/bench_areole.c $(AREOLE_PART_SRCS) include/ttypt/areole.h
	@mkdir -p bin
	${cc} ${CFLAGS} -o $@ tests/bench_areole.c src/areole.c ../libstele/src/stele_wq.c -lm -lpthread

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
