<script setup lang="ts">
// Mounts the interactive Chart.js figures of the current page. The page itself holds the
// static markup (<div data-leaf-app="...">); the application code, the vendored Chart.js
// and KaTeX, and the data written by docs/make.jl are served from the site's public files.
import { onBeforeUnmount, onMounted } from 'vue'
import { withBase } from 'vitepress'

declare global {
  interface Window { LeafCharts?: { mount(base: string): Promise<unknown>; unmount(): void } }
}

let loader: Promise<void> | undefined

function loadLoader(): Promise<void> {
  if (window.LeafCharts) return Promise.resolve()
  loader ??= new Promise((resolve, reject) => {
    const script = document.createElement('script')
    script.src = withBase('/leaf/loader.js')
    script.onload = () => resolve()
    script.onerror = () => reject(new Error(`could not load ${script.src}`))
    document.head.append(script)
  })
  return loader
}

onMounted(async () => {
  await loadLoader()
  await window.LeafCharts?.mount(withBase('/'))
})

onBeforeUnmount(() => window.LeafCharts?.unmount())
</script>

<template>
  <span class="leaf-charts-hook" aria-hidden="true"></span>
</template>
