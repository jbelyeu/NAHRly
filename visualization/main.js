// TODO: factor out code relating to chart into chart.js
// TODO: factor out code relating to table into new file called table.js

// You can only use import and export statements inside modules; not regular scripts.
// In some module systems, you can omit the file extension and the dot (e.g. '/modules/square').
// This doesn't work in native JavaScript modules.
import './chart.js'

// https://stackoverflow.com/a/39855320/6674256
import regions from '../data/sample_of_regions.js'

new Vue({
  el: '#app',
  vuetify: new Vuetify(),
  data () {
    return {
      regions: regions,
      names: Object.keys(regions),
      name: Object.keys(regions)[2],
      unique_colors: ['red', 'green', 'blue', 'magenta', 'cyan'],
      histogram_inferred_copy_number: 'all',
      bayes_inferred_copy_number: 'all',
      tippy: null,
      tippy_timer: null,
      wait: 1000
    }
  },
  // use regular functions, not arrow functions, when invoking "this" : https://michaelnthiessen.com/this-is-undefined/
  computed: {
    region () {
      return this.regions[this.name]
    },
    xs () {
      return this.region['samples'].map(sample => sample['normalized read depth'])
    },
    min_xs () {
      return Math.min(...this.xs)
    },
    max_xs () {
      return Math.max(...this.xs)
    },
    range_xs () {
      return this.max_xs - this.min_xs
    },
    possible_copy_numbers () {
      return this.region['possible copy numbers']
    },
    filter () {
      return (
        this.histogram_inferred_copy_number !== 'all'
        &&
        this.bayes_inferred_copy_number !== 'all'
      )
    },
    region_to_show () {
      if (this.filter) {
        // https://flaviocopes.com/how-to-clone-javascript-object/
        let region_to_show = _.cloneDeep(this.region)
        region_to_show['samples'] = this.samples(
          this.histogram_inferred_copy_number,
          this.bayes_inferred_copy_number
        )
        return region_to_show
      } else {
        return this.region
      }
    },
    // plotly only supports a limited number of html tags, e.g., "<br>"
    hover_texts_to_show () {
      return this.region_to_show['samples'].map(sample => {
        const first_line = 'histogram-inferred copy number: ' + sample['histogram-inferred copy number'] + '<br>'
        const second_line = 'possible copy numbers: ' + this.region_to_show['possible copy numbers'] + '<br>'
        const third_line = (
          'probabilities of possible copy numbers: ' +
          sample['probabilities of possible copy numbers'].map(value => value.toFixed(2)) +
          '<br>'
        )
        const fourth_line = 'bayes-inferred copy number: ' + sample['bayes-inferred copy number']
        return first_line + second_line + third_line + fourth_line
      })
    },
    xs_to_show () {
      return this.region_to_show['samples'].map(sample => sample['normalized read depth'])
    },
    ys_to_show () {
      return this.region_to_show['samples'].map(sample => sample['random number'])
    },
    colors_to_show () {
      return this.region_to_show['samples'].map(sample =>
        this.unique_colors[sample['histogram-inferred copy number'] - 1])
    },
    data () {
      return {
        traces: [
          {
            x: this.xs_to_show,
            y: this.ys_to_show,
            mode: 'markers',
            type: 'scatter',
            marker: {
              size: 12,
              color: this.colors_to_show
            },
            text: this.hover_texts_to_show,
            hoverinfo: 'text'
          }
        ],
        layout: {
          font: {
            family: 'Google Sans, sans-serif',
          },
          xaxis: {
            title: 'normalized read depth',
            range: [this.min_xs - 0.1*this.range_xs, this.max_xs + 0.1*this.range_xs]
          },
          // https://codepen.io/plotly/pen/KpLVzv ...
          yaxis: {
            showgrid: false,
            zeroline: false,
            showticklabels: false,
            range: [-0.2, 1.2]
          },
          hovermode: 'closest',
          hoverlabel: {
            bgcolor: '#fafafa'
          }
        }
      }
    }
  },
  methods: {
    samples (histogram_inferred_copy_number, bayes_inferred_copy_number) {
      return this.region['samples'].filter(sample => (
          sample['histogram-inferred copy number'] === histogram_inferred_copy_number
          &&
          sample['bayes-inferred copy number'] === bayes_inferred_copy_number
        )
      )
    },
    number_samples (histogram_inferred_copy_number, bayes_inferred_copy_number) {
      return this.samples(histogram_inferred_copy_number, bayes_inferred_copy_number).length
    },
    mouseenter (histogram_inferred_copy_number, bayes_inferred_copy_number) {
      this.histogram_inferred_copy_number = histogram_inferred_copy_number
      this.bayes_inferred_copy_number = bayes_inferred_copy_number
    },
    mouseleave () {
      this.histogram_inferred_copy_number = 'all'
      this.bayes_inferred_copy_number = 'all'
    },
    show_tippy () {
      this.tippy = tippy(this.$refs.table.$el, {
        content: 'Hover over cells of table to highlight samples in scatter plot below',
        placement: 'left',
        arrow: true,
        theme: 'google',
        trigger: 'manual',
        hideOnClick: true
      })
      this.tippy.show()
    }

  },
  mounted () {
    this.tippy_timer = setTimeout(this.show_tippy, this.wait)
  },
  beforeDestroy () {
    clearTimeout(this.tippy_timer)
    this.tippy.destroy(true)
  }

})
