require 'colorize'

def colorizedsh(cmd)
  puts cmd.yellow
  sh cmd
end

namespace :pkg do

  desc "Build doc"
  task :doc do
    puts "== Build doc".green
    output = `R --vanilla -e 'devtools::document(roclets=c("rd", "collate", "namespace"))'`
    output = output.gsub(/[Ee]rror/, "ERROR".red)
    output = output.gsub(/[Ww]arning/, "WARNING".yellow)
    puts output
  end

  desc "Install"
  task :install do
    puts "== Install".green
    sh 'R --vanilla -e "devtools::install()"'
  end

  desc "Test"
  task :test do
    puts "== Test".green
    sh 'R --vanilla -e "devtools::test()"'
  end

  desc "Run examples"
  task :run_examples do
    puts "== Run examples".green
    sh 'R --vanilla -e "devtools::run_examples()"'
  end

  desc "Check man"
  task :check_man do
    puts "== Check man".green
    sh 'R --vanilla -e "devtools::check_man()"'
  end

  desc "Check"
  task :check do
    puts "== Check".green
    sh 'R --vanilla -e "devtools::check()"'
  end

  namespace :site do

    vignettes =  FileList.new('vignettes/*.Rmd')
    desc "Build vignette for the website"
    task :build_vignette do
      puts "== Build vignette".green
      sh 'R --vanilla -e "pkgdown::build_articles()"'
    end

    desc "Build home of the website"
    task :build_home do
      puts "== Build home".green
      sh 'R --vanilla -e "pkgdown::build_home()"'
    end

    desc "Build reference of the website"
    task :build_reference do
      puts "== Build reference".green
      sh 'R --vanilla -e "pkgdown::build_reference()"'
    end

    desc "Make the pkg website"
    task :build do
      puts "== Build site".green
      sh 'R --vanilla -e "pkgdown::build_site()"'
    end

    desc "Open the pkg website"
    task :open do
      puts "== Open the site".green
      sh 'xdg-open docs/index.html'
    end

    desc "Remove docs dir"
    task :clean do
      puts "== Remove docs/ dir".green
      sh 'rm -rf docs/'
    end

    desc "Deploy the site on gh-pages"
    task :deploy => [:open] do
      puts "== Deploy to gh-pages".green
      colorizedsh "git clone .git gh-pages/"
      Dir.chdir("gh-pages/") do
        colorizedsh "git checkout gh-pages"
        colorizedsh 'rsync -ar ../docs/ .'
        colorizedsh "git add *"
        colorizedsh 'git commit -m "automatic deployment"'
        colorizedsh "git push origin gh-pages"
      end
      colorizedsh "rm -rf gh-pages/"
    end

  end

end
