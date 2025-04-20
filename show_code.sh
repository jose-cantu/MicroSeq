git ls-files \
    '*.py' '*.yml' '*.yaml' '*.toml' '*.md' '*.sh' \
    ':!:*.egg-info/*' ':!:*.ab1' ':!:LICENSES/*' \
| while read -r f; do
    echo -e "\n# =====  $f  =====\n"
    cat "$f"
done | less

